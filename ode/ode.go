// Package ode is a ordinary differential equation library that implements the
// Runge-Kutta methods https://en.wikipedia.org/wiki/Runge–Kutta_methods.
// The package requires a state space model, from the https://github.com/hammal/stateSpaceModel
// package, to describe the system.
package ode

import (
	"errors"
	"fmt"
	"math"
	"runtime"
	"sync"

	"github.com/hammal/adc/ssm"
	"gonum.org/v1/gonum/mat"
)

// DifferentiableSystem is an interface describing the necessary functions
// for a system to be solved using the ordinary differential equation solvers
// defined in this package.
type DifferentiableSystem interface {
	Derivative(t float64, state mat.Vector) mat.Vector
	Order() int
}

// RungeKutta holds the butcherTableau which describes the Runge Kutta method.
type RungeKutta struct {
	Description butcherTableau
}

// Compute ...
func (rk RungeKutta) Compute(from, to float64, value mat.Matrix, system DifferentiableSystem) (mat.Matrix, error) {
	M, N := value.Dims()

	res := make([]chan mat.Vector, N)

	var wg sync.WaitGroup

	wg.Add(N)

	tmpVec := value.(*mat.Dense)
	resValue := mat.NewDense(M, N, nil)

	for column := 0; column < N; column++ {
		res[column] = make(chan mat.Vector)

		go func(
			from, to float64,
			initalValue mat.Vector,
			system DifferentiableSystem,
			wg *sync.WaitGroup,
			receiveChannel chan mat.Vector,
		) {
			// Report when done
			defer wg.Done()
			// Close receive channel when done
			defer close(receiveChannel)
			// Compute the adaptive vector function
			res, _ := rk.computeVec(from, to, initalValue, system)
			// Send back the result through channels
			receiveChannel <- res
		}(from, to, tmpVec.ColView(column), system, &wg, res[column])
	}

	var tmp mat.Vector

	for colindex, col := range res {
		tmp = <-col
		for rowIndex := 0; rowIndex < M; rowIndex++ {
			resValue.Set(rowIndex, colindex, tmp.AtVec(rowIndex))
		}
	}

	// This might not serve any purpose
	// wg.Wait()

	return resValue, nil
}

// computeVec computes the update for a Runge-Kutta system based on a current value at t = from
// , a target time t = to, a initial value x(t=from) = value and a system model ode.
// The function returns a result and associated error vector
func (rk RungeKutta) computeVec(from, to float64, value mat.Vector, system DifferentiableSystem) (mat.Vector, mat.Vector) {

	// Define variables
	var (
		tempV *mat.VecDense
	)

	// State order
	M, _ := value.Dims()

	// The precomputed derivative points
	K := make([]mat.Vector, rk.Description.stages)

	// Step length
	h := to - from

	// Compute all derivative points
	for index := range K {

		// Check system type for different implementations
		switch sys := system.(type) {

		// Linear state space models:
		// The state can be computed in closed form such x(t)= e^(A (to- from)) x(0)
		// And thus only the inputs need to be solved using the Runge Kutta methods.
		// Note that the state is added in a later segment
		case *ssm.LinearStateSpaceModel:
			tmpM := sys.StateSpaceOrder()
			tempV = mat.NewVecDense(tmpM, nil)

		default:
			// Initialize an intermediate vector
			// fmt.Printf("This should not happen for %T", system)
			tempV = mat.NewVecDense(M, nil)
			tempV.CloneVec(value)
		}
		// Compute the relevant vector by combining previously computed derivate points
		// according to Butcher Tableau.
		for index2, a := range rk.Description.rungeKuttaMatrix[index] {
			tempV.AddScaledVec(tempV, h*a, K[index2])
		}
		// Insert the new derivate point
		// These can be implemented differently depending on underlying model
		K[index] = system.Derivative(from+h*rk.Description.nodes[index], tempV)
	}

	// Summarise the different derivation points
	switch sys := system.(type) {
	// For the Linear state space model this can be done in closed form
	case *ssm.LinearStateSpaceModel:
		// compute e^(AT_s) and move state forward
		var tmpMatrix mat.Dense
		tmpMatrix.Scale(to-from, sys.A)
		// fmt.Printf("A T_s = \n%v\n", mat.Formatted(&tmpMatrix))
		// TEMP This solution is quite sensitive when the condition number of the
		// matrix tmpMatrix above is large.
		tmpMatrix.Exp(&tmpMatrix)

		// Initialize tempV as the initial state with the state dynamics applied.
		tempV.Reset()
		tempV.MulVec(&tmpMatrix, value)
		// fmt.Printf("e^AT_s X(0) = \n%v\n", mat.Formatted(tempV))
	default:
		// Reset tempV to the initial value
		tempV.CloneVec(value)
	}

	// Initialize the error vector
	err := mat.NewVecDense(M, nil)
	// Sum up the different contributions with relevant weights.
	for index, k := range K {
		tempV.AddScaledVec(tempV, h*rk.Description.weights[0][index], k)
		// If the Butcher Tableau allows for adaptive error computation
		if len(rk.Description.weights) == 2 {
			err.AddScaledVec(err, h*(rk.Description.weights[1][index]-rk.Description.weights[0][index]), k)
		}
	}

	// fmt.Printf("Result from computeVec\n%v\n", value)

	return tempV, err
}

// AdaptiveCompute implements an adaptive version which for a
// given error tolerance err. Makes recursive steps such that the local error
// never exceeds the error specification. The functions returns a result matrix
// and a error struct.
func (rk RungeKutta) AdaptiveCompute(from, to, errorTarget float64, value mat.Matrix, system DifferentiableSystem) (mat.Matrix, error) {
	M, N := value.Dims()

	res := make([]chan mat.Vector, N)
	err := make([]chan error, N)

	var wg sync.WaitGroup

	wg.Add(N)

	tmpVec := value.(*mat.Dense)
	resValue := mat.NewDense(M, N, nil)

	for column := 0; column < N; column++ {
		res[column] = make(chan mat.Vector, 1)
		err[column] = make(chan error, 1)

		go func(
			from, to, errorTarget float64,
			initalValue mat.Vector,
			system DifferentiableSystem,
			wg *sync.WaitGroup,
			receiveChannel chan mat.Vector,
			errorChannel chan error,
		) {
			// Report when done
			defer wg.Done()
			// Close receive channel when done
			defer close(receiveChannel)
			defer close(errorChannel)
			// Compute the adaptive vector function
			res, err := rk.adaptiveComputeVec(from, to, errorTarget, initalValue, system)
			// Send back the result through channels
			receiveChannel <- res
			errorChannel <- err
		}(from, to, errorTarget, tmpVec.ColView(column), system, &wg, res[column], err[column])
	}

	// Check for errors and return immediately if so.
	for _, e := range err {
		if ee := <-e; ee != nil {
			return nil, ee
		}
	}

	var tmp mat.Vector

	for colindex, col := range res {
		tmp = <-col
		for rowIndex := 0; rowIndex < M; rowIndex++ {
			resValue.Set(rowIndex, colindex, tmp.AtVec(rowIndex))
		}
	}

	// wg.Wait()

	return resValue, nil
}

// adaptiveComputeVec computes...
func (rk RungeKutta) adaptiveComputeVec(from, to, err float64, value mat.Vector, system DifferentiableSystem) (mat.Vector, error) {
	var (
		tmpVec                   mat.Vector
		tmpState1                *mat.VecDense
		tmpState2                []*mat.VecDense
		currentError             float64
		tnow                     float64
		tnext, stepLength        []float64
		count                    int
		NumberOfParallelSearches int
		WinningCandidate         int
	)
	// Set max number of iterations
	const maxNumberOfIterations int = 1000
	const searchStepFactor = 2

	// Parallel Receive channels
	NumberOfParallelSearches = 1
	if runtime.NumCPU() > 1 {
		NumberOfParallelSearches = runtime.NumCPU() * 8
	}
	ResponseChannel := make([]chan mat.Vector, NumberOfParallelSearches)
	ComChannel := make([]chan int, NumberOfParallelSearches)

	for channelIndex := range ResponseChannel {
		ResponseChannel[channelIndex] = make(chan mat.Vector)
		ComChannel[channelIndex] = make(chan int)
	}

	tnext = make([]float64, NumberOfParallelSearches)
	stepLength = make([]float64, NumberOfParallelSearches)

	// fmt.Printf("Number of processors: %v\n", NumberOfParallelSearches)

	// Initialize current time
	tnow = from
	for index := range tnext {
		tnext[index] = tnow
		// Set the update rate
		stepLength[index] = to - from
	}

	M := value.Len()
	tmpState1 = mat.NewVecDense(M, nil)
	tmpState1.CloneVec(value)
	tmpState2 = make([]*mat.VecDense, NumberOfParallelSearches)
	for index := range tmpState2 {
		tmpState2[index] = mat.NewVecDense(M, nil)
	}

	// Repeat until time to is reached
	// fmt.Println()
	for tnow < to {
		// Set target time
		// fmt.Printf("\nNew Round! Tnext: ")
		for index := range tnext {
			tnext[index] = math.Min(to, tnow+stepLength[index])
			// fmt.Printf("%v, ", tnext[index])
		}
		// fmt.Println()
		// fmt.Printf("tnow: %v\n", tnow)
		// Repeat until target error is reached
		for true {
			// Copy the current state into tmpState
			for index := range tnext {
				tmpState2[index].CopyVec(tmpState1)
				// Execute the Runge Kutta computation
				// fmt.Printf("Started index %v\n", index)
				go func(respChannel chan mat.Vector, com chan int, tnow, tnext float64, currentState *mat.VecDense) {
					tmpVec, currentErrorVector := rk.computeVec(tnow, tnext, currentState, system)
					for {
						select {
						case respChannel <- currentErrorVector:
							respChannel <- tmpVec
							return
						case <-com:
							return
						}
					}

				}(ResponseChannel[index], ComChannel[index], tnow, tnext[index], tmpState2[index])
			}

			// Find the WinningCandidate. Longest stepLength take presidence
			WinningCandidate = -1
			for index := range tnext {
				// fmt.Printf("Waiting for index %v\n", index)
				errVec := <-ResponseChannel[index]
				tmpVec = <-ResponseChannel[index]
				currentError := mat.Norm(errVec, 1)
				// fmt.Printf("Received errorVector \n%v\n", mat.Formatted(errVec))
				if currentError < err {
					WinningCandidate = index
					// close go routines
					for index2 := WinningCandidate + 1; index2 < len(tnext); index2++ {
						ComChannel[index2] <- 1
					}
					break
				}
			}

			if WinningCandidate >= 0 {
				break
			}

			// If not
			// Increment counter and check if we are allowed more trials
			count++
			if count >= maxNumberOfIterations {
				errorString := fmt.Sprintf("Maximum number of iterations reached adaptive Runge-Kutta doesn't converge\n last error was = %v, and only %2.f percent of time was computed", currentError, tnext[0]/to*100)
				return nil, errors.New(errorString)
			}
			// Half the next integration interval and try again
			// fmt.Print("Reduced steplength: ")
			for index := range stepLength {
				stepLength[index] = (tnext[NumberOfParallelSearches-1] - tnow) * math.Pow(searchStepFactor, -float64(index+1))
				tnext[index] = stepLength[index] + tnow
				// fmt.Printf("%v, ", stepLength[index])
			}
			// fmt.Println()

		}
		// Save this state and update tnow
		tmpState1.CopyVec(tmpVec)
		// fmt.Printf("%v, ", WinningCandidate)
		// fmt.Printf("Winning Candidate %v\n", WinningCandidate)
		tnow = tnext[WinningCandidate]
		tmpStepLength := stepLength[WinningCandidate]
		// fmt.Printf("time: %.10e\n", tnow)
		// fmt.Print("New StepLengths: ")
		for index := range stepLength {
			stepLength[index] = tmpStepLength * math.Pow(searchStepFactor, float64(NumberOfParallelSearches/2-index))
			// fmt.Printf("%v, ", stepLength[index])
		}
		count = 0
		// fmt.Println()

	}

	// Successful integration!  Return nil error
	return tmpState1, nil
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
		{-8. / 27., 2., -3544. / 2565., 1859. / 4104., -11. / 40.},
	}
	rk := RungeKutta{temp}
	return &rk
}
