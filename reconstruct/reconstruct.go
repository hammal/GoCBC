package reconstruct

type Reconstruction interface {
	// Runs the reconstruction and returns the result
	Reconstruct() [][]float64
}
