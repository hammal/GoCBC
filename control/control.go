package control

import (
	"gonum.org/v1/gonum/mat"
)

// Control interface holds the Control instance which is capable of simulating
// the system as well as providing the filter contributions of the controls at a given index.
type Control interface {
	// Executes the control simulation and return the observations
	Simulate() [][]float64
	// Get the control contribution for simulation at index
	getControlSimulationContribution(index int) (mat.Vector, error)
	// Get the control contribution for filtering at index
	GetForwardControlFilterContribution(index int) (mat.Vector, error)
	GetBackwardControlFilterContribution(index int) (mat.Vector, error)
	// Precompute filter decisions based on filter dynamics
	PreComputeFilterContributions(forwardDynamics, backwardDynamics mat.Matrix)
	// Length of control
	GetLength() int
	// get the sample period
	GetTs() float64
}

// Cache interface is an abstraction that is heavily used when precomputing
// lookup tables.
type ControlVector interface {
	GetVector(uint) mat.Vector
}
