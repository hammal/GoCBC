// Package ode is a ordinary differential equation library that implements the
// Runge-Kutta methods https://en.wikipedia.org/wiki/Runge–Kutta_methods.
// The package requires a state space model, from the https://github.com/hammal/stateSpaceModel
// package, to describe the system.
package ode

import (
	"errors"
	"math"
	"sync"

	"github.com/hammal/adc/ssm"
	"gonum.org/v1/gonum/mat"
)

type DifferentiableSystem interface {
	Derivative(t float64, state mat.Vector) mat.Vector
}

// RungeKutta holds the butcherTableau which describes the Runge Kutta method.
type RungeKutta struct {
	Description butcherTableau
}

func (rk RungeKutta) Compute(from, to float64, value mat.Matrix, system DifferentiableSystem) mat.Matrix {
	M, N := value.Dims()

	res := make([]mat.Vector, N)
	err := make([]mat.Vector, N)

	var wg sync.WaitGroup

	wg.Add(N)

	for column := 0; column < N; column++ {
		// res[column] = mat.NewVecDense(M, nil)
		v := value.(*mat.Dense)
		res[column] = v.ColView(column)

		go func(column int, err []mat.Vector, from, to float64, res []mat.Vector, system DifferentiableSystem, wg *sync.WaitGroup) {
			err[column] = rk.computeVec(from, to, res[column], system, wg)
		}(column, err, from, to, res, system, &wg)
	}

	wg.Wait()

	errMatrix := mat.NewDense(M, N, nil)
	for row := 0; row < M; row++ {
		for column := 0; column < N; column++ {
			v := value.(*mat.Dense)
			v.Set(row, column, res[column].AtVec(row))
			errMatrix.Set(row, column, err[column].AtVec(row))
		}
	}
	return errMatrix
}

// computeVec the update for a Runge-Kutta system based on a current value at t = from
// , a target time t = to, a initial value x(t=from) = value and a system model ode.
// When algorithm is finished the result is copied into the value.
func (rk RungeKutta) computeVec(from, to float64, value mat.Vector, system DifferentiableSystem, sync *sync.WaitGroup) mat.Vector {
	defer sync.Done()
	// Variables
	var (
		tempV mat.VecDense
	)

	// State order
	M, _ := value.Dims()
	// The precomputed derivative points
	K := make([]mat.Vector, rk.Description.stages)
	// Step length
	h := to - from
	for index := range K {
		switch sys := system.(type) {
		case *ssm.LinearStateSpaceModel:
			tempV = *mat.NewVecDense(M, nil)
			for _, inp := range sys.Input {
				tempV.AddVec(&tempV, inp.Value(from+h*rk.Description.nodes[index]))
			}

			K[index] = &tempV
		default:
			// Initialize an intermediate vector
			// tempV := mat.NewVecDense(M, nil)
			tempV.CloneVec(value)
			// Compute the relevant vector by combining previously computed derivate points
			// according to Butcher Tableau.
			for index2, a := range rk.Description.rungeKuttaMatrix[index] {
				tempV.AddScaledVec(&tempV, h*a, K[index2])
			}
			// Insert the new derivate point
			// These can be implemented differently depending on underlying model
			K[index] = system.Derivative(from+h*rk.Description.nodes[index], &tempV)
		}

	}
	switch sys := system.(type) {
	case *ssm.LinearStateSpaceModel:
		// compute e^(AT_s) and move state forward
		var tmpMatrix mat.Dense
		tmpMatrix.Scale(to-from, sys.A)
		tmpMatrix.Exp(&tmpMatrix)
		tempV.MulVec(&tmpMatrix, value)
	default:
		// Reset tempV
		tempV.CloneVec(value)
	}

	// Initialize the error vector
	err := mat.NewVecDense(M, nil)
	// Sum up the different contributions with relevant weights.
	for index, k := range K {
		tempV.AddScaledVec(&tempV, h*rk.Description.weights[0][index], k)
		// If the Butcher Tableau allows for adaptive error computation
		if len(rk.Description.weights) == 2 {
			err.AddScaledVec(err, h*(rk.Description.weights[1][index]-rk.Description.weights[0][index]), k)
		}
	}

	value = &tempV
	return err
}

func (rk RungeKutta) AdaptiveCompute(from, to, errorTarget float64, value mat.Matrix, system DifferentiableSystem) error {
	M, N := value.Dims()

	res := make([]mat.Vector, N)
	err := make([]error, N)

	var wg sync.WaitGroup

	wg.Add(N)

	for column := 0; column < N; column++ {
		// res[column] = mat.NewVecDense(M, nil)
		v := value.(*mat.Dense)
		res[column] = v.ColView(column)

		go func(column int, err []error, from, to, errorTarget float64, res []mat.Vector, system DifferentiableSystem, wg *sync.WaitGroup) {
			err[column] = rk.adaptiveComputeVec(from, to, errorTarget, res[column], system, wg)
		}(column, err, from, to, errorTarget, res, system, &wg)
	}

	wg.Wait()

	for column := 0; column < N; column++ {
		for row := 0; row < M; row++ {
			v := value.(*mat.Dense)
			v.Set(row, column, res[column].AtVec(row))
		}
		if err[column] != nil {
			return err[column]
		}
	}
	return nil
}

// AdaptiveCompute implements an adaptive version which for a
// given error tolerance err. Makes recursive steps such that the local error
// never exceeds the error specification.
func (rk RungeKutta) adaptiveComputeVec(from, to, err float64, value mat.Vector, system DifferentiableSystem, thisSync *sync.WaitGroup) error {
	defer thisSync.Done()
	var (
		tmpState1          *mat.VecDense
		tmpState2          *mat.VecDense
		currentErrorVector mat.Vector
		currentError       float64
		tnow, tnext        float64
		count              int
	)
	// Set max number of iterations
	const maxNumberOfIterations int = 10000

	// Initialize current time
	tnow = from

	M := value.Len()
	tmpState1 = mat.NewVecDense(M, nil)
	tmpState2 = mat.NewVecDense(M, nil)
	tmpState1.CloneVec(value)

	// Repeat until time to is reached
	for tnow < to {
		// Set target time
		tnext = to
		// Repeat until target error is reached
		for true {
			// Copy the current state into tmpState
			tmpState2.CopyVec(tmpState1)
			// Execute the Runge Kutta computation
			var sync sync.WaitGroup
			sync.Add(1)
			currentErrorVector = rk.computeVec(tnow, tnext, tmpState2, system, &sync)
			// Reset and compute Error
			currentError = 0.
			// count = 0
			for index := 0; index < tmpState2.Len(); index++ {
				currentError += math.Abs(currentErrorVector.AtVec(index))
			}
			// Has the target error been achived?
			// fmt.Printf("Err %v, time %v\n", currentError, tnow)
			if currentError < err {
				break
			}
			// Half the next integration interval and try again
			tnext = (tnext-tnow)/2. + tnow

			// Increment counter and check if we are allowed more trials
			count++
			if count >= maxNumberOfIterations {
				return errors.New("Maximum number of iterations reached adaptive Runge-Kutta doesn't converge")
			}
		}
		// Save this state and update tnow
		tmpState1.CopyVec(tmpState2)
		tnow = tnext

	}
	value = tmpState1

	// Successful integration!  Return nil error
	return nil
}

// TODO: Implement Compute and Adaptive Compute for Matrices of states

// NewRK4 function returns a forth order Runge-Kutta object
func NewRK4() *RungeKutta {
	var temp butcherTableau
	temp.stages = 4
	temp.nodes = []float64{0, 1. / 2., 1. / 2., 1}
	temp.weights = [][]float64{{1. / 6., 1. / 3., 1. / 3., 1. / 6.}}
	temp.rungeKuttaMatrix = [][]float64{
		nil,
		{1. / 2.},
		{0, 1. / 2.},
		{0, 0, 1.},
	}
	rk := RungeKutta{temp}
	return &rk
}

// NewEulerMethod returns a pointer to a Runge-Kutta that does the Euler method.
func NewEulerMethod() *RungeKutta {
	var temp butcherTableau
	temp.stages = 1
	temp.nodes = []float64{0}
	temp.weights = [][]float64{{1}}
	rk := RungeKutta{temp}
	return &rk
}

// butcherTableau which describes the approximate solution, see https://en.wikipedia.org/wiki/Runge–Kutta_methods.
type butcherTableau struct {
	stages           int
	weights          [][]float64
	nodes            []float64
	rungeKuttaMatrix [][]float64
}

// NewFehlberg45 implements https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method
func NewFehlberg45() *RungeKutta {
	var temp butcherTableau
	temp.stages = 6
	temp.nodes = []float64{0, 1. / 4., 3. / 8., 12. / 13., 1., 1. / 2.}
	temp.weights = [][]float64{
		{16. / 135., 0, 6656. / 12825., 28561. / 56430., -9. / 50., 2. / 55.},
		{25. / 216., 0, 1408. / 2565., 2197. / 4104., -1. / 5., 0},
	}
	temp.rungeKuttaMatrix = [][]float64{
		nil,
		{1. / 4.},
		{3. / 32., 9. / 32.},
		{1932. / 2197., -7200. / 2197., 7296. / 2197.},
		{439. / 216., -8., 3680. / 513., -845. / 4104.},
		{-8. / 27., 2, -3544. / 2565., 1859. / 4104., -11. / 40.},
	}
	rk := RungeKutta{temp}
	return &rk
}
