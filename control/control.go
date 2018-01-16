package control

import (
	"github.com/hammal/adc/ode"
	"github.com/hammal/adc/signal"
	"github.com/hammal/adc/ssm"
	"gonum.org/v1/gonum/mat"
)

// Control interface holds the Control instance which is capeable of simulating
// the system as well as providing the netto controls at a given index.
type Control interface {
	// Executes the control simulation and return the observations
	Simulate() [][]float64
	// Get the control contribution for simulation at index
	getControlSimulationContribution(index int) []mat.Vector
	// Get the control contribution for filtering at index
	GetForwardControlFilterContribution(index int) []mat.Vector
	GetBackwardControlFilterContribution(index int) []mat.Vector
	// Precompute filter descisions
	PreComputeFilterContributions(eta2 []float64, Vf mat.Vector, Vb mat.Vector)
	// Length of control
	Length() int
}

// zeroOrderHold solves the problem
// int_0^T_s e^(A(T_s - t)) dt
// by solving the initial value problem
// x'(t) = A x(t) + B
func zeroOrderHold(A *mat.Dense, t float64) mat.Matrix {
	// Determine allowed error
	const err float64 = 1e-9

	// Define variables
	M, _ := A.Dims()
	res := mat.NewDense(M, M, nil)
	input := make([]signal.VectorFunction, 1)

	// For each possible unity input
	for column := 0; column < M; column++ {
		// Create unit vector
		tmp := mat.NewVecDense(M, nil)
		tmp.SetVec(column, 1.)
		input[0] = signal.NewInput(func(arg1 float64) float64 { return 1. }, tmp)
		// Assign temporary state space model and initial state vector
		tempSSM := ssm.NewLinearStateSpaceModel(A, A, input)
		tempState := mat.NewVecDense(M, nil)
		// Solve initial value problem using adaptive Runge-Kutta method.
		od := ode.NewFehlberg45()
		od.AdaptiveCompute(0., t, err, tempState, tempSSM)

		// Fill up result matrix for given unit vector.
		for row := 0; row < M; row++ {
			res.Set(row, column, tempState.AtVec(row))
		}
	}
	return res
}
