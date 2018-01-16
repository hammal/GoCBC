package ssm

import (
	"github.com/hammal/adc/signal"
	"gonum.org/v1/gonum/mat"
)

type MemoryLessLinearStateSpaceModel struct {
	// Observation matrix
	C *mat.Dense
	// input
	Input signal.VectorFunction
}

// StateDerivative returns the state derivative.
// x'(t) = Ax(t) + Bu(t)
// where state = x(t) at an arbitrary time t. Furthermore, Bu is the input vector field.
func (model MemoryLessLinearStateSpaceModel) StateDerivative(t float64, state *mat.VecDense) *mat.VecDense {
	return model.Input.Bu(t)
}

// StateObservation returns the observed stateDerivate
// y(t) = C 0
// where
// state = x(t) and t is an arbitrary time.
func (model MemoryLessLinearStateSpaceModel) StateObservation(t float64, state *mat.VecDense) *mat.VecDense {
	n, _ := model.C.Dims()
	return mat.NewVecDense(n, nil)
}

func (ssm MemoryLessLinearStateSpaceModel) StateSpaceOrder() int {
	m, _ := ssm.Input.B.Dims()
	return m
}

func (ssm MemoryLessLinearStateSpaceModel) ObservationSpaceOrder() int {
	m, _ := ssm.C.Dims()
	return m
}
