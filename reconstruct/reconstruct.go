package reconstruct

import "gonum.org/v1/gonum/mat"

type Reconstruction interface {
	// Runs the reconstruction and returns the result
	Reconstruct() [][]*float64

	// Get State dynamics, used to pre-compute control decision vectors.
	GetStateDynamics() (Af, Ab *mat.Dense)
}
