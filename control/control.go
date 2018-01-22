package control

import (
	"fmt"

	"github.com/hammal/adc/ode"
	"github.com/hammal/adc/signal"
	"github.com/hammal/adc/ssm"
	"gonum.org/v1/gonum/mat"
)

// Control interface holds the Control instance which is capable of simulating
// the system as well as providing the filter contributions of the controls at a given index.
type Control interface {
	// Executes the control simulation and return the observations
	Simulate() [][]float64
	// Get the control contribution for simulation at index
	getControlSimulationContribution(index int) []mat.Vector
	// Get the control contribution for filtering at index
	GetForwardControlFilterContribution(index int) []mat.Vector
	GetBackwardControlFilterContribution(index int) []mat.Vector
	// Precompute filter decisions based on filter dynamics
	PreComputeFilterContributions(forwardDynamics, backwardDynamics mat.Matrix)
	// Length of control
	GetLength() int
	// get the sample period
	GetTs() float64
}

// zeroOrderHold solves the problem
// int_0^T_s e^(A(T_s - t)) dt
// by solving the initial value problem
// x'(t) = A x(t) + B
func zeroOrderHold(A mat.Matrix, t float64) mat.Matrix {
	// Determine allowed error
	const err float64 = 1e-9

	fmt.Println(mat.Formatted(A))

	// Define variables
	M, _ := A.Dims()
	res := mat.NewDense(M, M, nil)
	input := make([]signal.VectorFunction, 1)
	od := ode.NewFehlberg45()

	// For each possible unity input
	for column := 0; column < M; column++ {
		// Create unit vector
		tmp := mat.NewVecDense(M, nil)
		tmp.SetVec(column, 1.)
		input[0] = signal.NewInput(func(arg1 float64) float64 { return 1. }, tmp)
		// Assign temporary state space model and initial state vector
		tempSSM := ssm.NewLinearStateSpaceModel(A, A, input)
		tempState := mat.NewDense(M, 1, nil)
		// Solve initial value problem using adaptive Runge-Kutta Fehlberg4(5) method.
		// od.AdaptiveCompute(0., t, err, tempState, tempSSM)
		tempRes, _ := od.Compute(0., t, tempState, tempSSM)

		// fmt.Println(mat.Formatted(tempRes))
		// fmt.Println(M)

		// Fill up result matrix for given unit vector.
		for row := 0; row < M; row++ {
			res.Set(row, column, tempRes.At(row, 0))
		}
	}
	return res
}
