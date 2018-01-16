// Package ode is a ordinary differential equation library that implements the
// Runge-Kutta methods https://en.wikipedia.org/wiki/Runge–Kutta_methods.
// The package requires a state space model, from the https://github.com/hammal/stateSpaceModel
// package, to describe the system.
package ode

import (
	"errors"
	"math"

	"github.com/hammal/adc/ssm"
	"gonum.org/v1/gonum/mat"
)

// RungeKutta holds the butcherTableau which describes the Runge Kutta method.
type RungeKutta struct {
	Description butcherTableau
}

// Compute the update for a Runge-Kutta system based on a current value at t = from
// , a target time t = to, a initial value x(t=from) = value and a system model ode.
// When algorithm is finished the result is copied into the value.
func (rk RungeKutta) Compute(from, to float64, value *mat.VecDense, ode ssm.StateSpaceModel) *mat.VecDense {
	// State order
	M, _ := value.Dims()
	// The precomputed derivative points
	K := make([]*mat.VecDense, rk.Description.stages)
	// Step length
	h := to - from
	for index := range K {
		// Initialize an intermediate vector
		tempV := mat.NewVecDense(M, nil)
		tempV.CopyVec(value)
		// Compute the relevant vector by combining previously computed derivate points
		// according to Butcher Tableau.
		for index2, a := range rk.Description.rungeKuttaMatrix[index] {
			tempV.AddScaledVec(tempV, h*a, K[index2])
		}
		// Insert the new derivate point
		K[index] = ode.StateDerivative(from+h*rk.Description.nodes[index], tempV)
	}
	// Initialize the error vector
	err := mat.NewVecDense(M, nil)
	// Sum up the different contributions with relevant weights.
	for index, k := range K {
		value.AddScaledVec(value, h*rk.Description.weights[0][index], k)
		// If the Butcher Tableau allows for adaptive error computation
		if len(rk.Description.weights) == 2 {
			err.AddScaledVec(err, h*(rk.Description.weights[1][index]-rk.Description.weights[0][index]), k)
		}
	}
	return err
}

// This adaptive Runge Kutta version implies an
func (rk RungeKutta) AdaptiveCompute(from, to, err float64, value *mat.VecDense, ode ssm.StateSpaceModel) error {
	var (
		tmpState           *mat.VecDense
		currentErrorVector *mat.VecDense
		currentError       float64
		tnow, tnext        float64
		count              int
	)
	// Set max number of iterations
	const maxNumberOfIterations int = 10000

	// Initialize current time
	tnow = from

	M, _ := value.Dims()
	tmpState = mat.NewVecDense(M, nil)

	// Repeat until time to is reached
	for tnow < to {
		// Set target time
		tnext = to
		// Repeat until target error is reached
		for true {
			// Copy the current state into tmpState
			tmpState.CopyVec(value)
			// Execute the Runge Kutta computation
			currentErrorVector = rk.Compute(tnow, tnext, tmpState, ode)
			// Reset and compute Error
			currentError = 0.
			// count = 0
			for index := 0; index < tmpState.Len(); index++ {
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
			count += 1
			if count >= maxNumberOfIterations {
				return errors.New("Max number of iterations reached adaptive Runge-Kutta doesn't converge.")
			}
		}
		// Save this state and update tnow
		value.CopyVec(tmpState)
		tnow = tnext

	}
	// Successful integration!  Return nil error
	return nil
}

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