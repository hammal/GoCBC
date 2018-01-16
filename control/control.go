package control

import "gonum.org/v1/gonum/mat"

// Control interface holds the Control instance which is capeable of simulating
// the system as well as providing the netto controls at a given index.
type Control interface {
	// Executes the control simulation and return the observations
	Simulate() [][]float64
	// Get the control contribution for simulation at index
	getControlSimulationContribution(index int) []mat.Vector
	// Get the control contribution for filtering at index
	GetCotnrolFilterContribution(index int) []mat.Vector
	// Precompute filter descisions
	PreComputeFilterContributions(eta2 []float64, Vf mat.Vector, Vb mat.Vector)
}

func zeroOrderHold(m *mat.Dense, t float64) mat.Matrix {
	m.Scale(t, m)
	var tmp mat.Dense
	tmp.Exp(m)
	return &tmp
}
