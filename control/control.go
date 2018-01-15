package control

import "gonum.org/v1/gonum/mat"

// Control interface holds the Control instance which is capeable of simulating
// the system as well as providing the netto controls at a given index.
type Control interface {
	// Executes the control simulation and return the observations
	Simulate() [][]float64
	// Get the control contribution at index
	GetControlContribution(index int) []mat.VecDense
}
