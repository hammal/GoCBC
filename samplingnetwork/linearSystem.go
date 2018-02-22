package samplingnetwork

import (
	"github.com/hammal/adc/signal"
	"github.com/hammal/adc/ssm"
	"gonum.org/v1/gonum/mat"
)

type LinearSystem struct {
	A mat.Matrix
	B mat.Matrix
	C mat.Matrix
}

func (sys LinearSystem) InputSpaceOrder() int {
	_, N := sys.B.Dims()
	return N
}

func (sys LinearSystem) OutputSpaceOrder() int {
	M, _ := sys.C.Dims()
	return M
}

func (sys LinearSystem) StateSpaceOrder() int {
	M, _ := sys.A.Dims()
	return M
}

// LinearSystemToLinearStateSpaceModel converts a system and an array of input
// functions into a linear state space model.
func LinearSystemToLinearStateSpaceModel(system LinearSystem, inputFunction []func(float64) float64) *ssm.LinearStateSpaceModel {

	if len(inputFunction) != system.InputSpaceOrder() {
		panic("The B vector and size of inputFunction don't agree")
	}

	input := make([]signal.VectorFunction, system.InputSpaceOrder())

	var tmpB mat.Dense

	tmpB.Clone(system.B)

	for index := range inputFunction {
		input[index] = signal.NewInput(inputFunction[index], tmpB.ColView(index))
	}

	return ssm.NewLinearStateSpaceModel(system.A, system.C, input)
}
