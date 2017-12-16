package simulator

import "gonum.org/v1/gonum/mat"

type Input interface {
	// At is a function can be evaluated at time t for a given argument.
	At(t float64) *mat.Vector
}

type Signal struct {
	signal func(float64) float64
	vector *mat.VecDense
	kind   string
}

// The controls are just signals
type Control Signal
