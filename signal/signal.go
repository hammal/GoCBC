package signal

import (
	"math"

	"gonum.org/v1/gonum/mat"
)

// Signal holds the signal interface
type Signal interface {
	Value(float64) mat.Vector
}

// DiracDelta is a Dirac delta distribution as defined in
// https://en.wikipedia.org/wiki/Dirac_delta_function
func DiracDelta(x float64) float64 {
	// These could all be done offline
	var a = 1e-9
	var a2 = a * a
	var C1 = 1. / (math.Abs(a) * math.Sqrt(math.Pi))
	// Return distribution value at x
	return C1 * math.Exp(-x*x/a2)
}
