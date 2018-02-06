package ssm

import (
	"gonum.org/v1/gonum/mat"
)

// StateSpaceModel interface has three parts:
//
// 1) The f function which returns the differential state evaluated at time t
// and state(t).
//
// 2) The system impulse response h(t) evaluated at N points t(0)...t(N-1).
//
// 3) Similarly the frequencyResponse H(f) evaluated at f(0)...f(N-1)
type StateSpaceModel interface {
	// This is the derivative of a state space mordinaryDifferentialEquationl
	Derivative(t float64, state mat.Vector) mat.Vector
	// This is the observedState
	Observation(t float64, state mat.Vector) mat.Vector
	// Returns the state space order
	StateSpaceOrder() int
	// Returns the observation space order.
	ObservationSpaceOrder() int
	// Returns the input space Order
	InputSpaceOrder() int
	// General Order
	Order() int
}
