package ssm

import (
	"gonum.org/v1/gonum/mat"
)

// SSM interface has three parts:
//
// 1) The f function which returns the differential state evaluated at time t
// and state(t).
//
// 2) The system impulse response h(t) evaluated at N points t(0)...t(N-1).
//
// 3) Similarly the frequencyResponse H(f) evaluated at f(0)...f(N-1)
type SSM interface {
	// This is the derivative of a state space model
	StateDerivative(t float64, state *mat.VecDense) *mat.VecDense
	// This is the observedState
	StateObservation(t float64, state *mat.VecDense) *mat.VecDense
	// Impulse response for system
	// ImpulseResponse(t []float64) [][][]float64
	// Frequency response
	// FrequencyResponse(f []float64) [][][]complex128
}
