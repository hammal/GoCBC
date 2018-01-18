package reconstruct

import (
	"fmt"
	"testing"

	"github.com/hammal/adc/signal"
	"github.com/hammal/adc/ssm"
	"gonum.org/v1/gonum/mat"
)

func TestCare(t *testing.T) {
	N := 4
	stageGain := 10.
	sigma_u := 1e-5
	sigma_z := 1e-3
	B := mat.NewVecDense(N, []float64{1, 0, 0, 0})
	input := make([]signal.VectorFunction, 1)
	input[0] = signal.NewInput(func(arg1 float64) float64 {
		return 1.
	}, B)
	stateSpaceModel := ssm.NewIntegratorChain(N, stageGain, input)
	Q := mat.NewDense(N, N, nil)
	Q.Outer(sigma_u, B, B)

	data := make([]float64, 1)
	data[0] = 1. / sigma_z
	Rinv := mat.NewDense(1, 1, data)

	Vf := mat.NewDense(N, N, nil)
	for index := 0; index < N; index++ {
		Vf.Set(index, index, 1.)
	}
	care(stateSpaceModel.A, stateSpaceModel.C, Rinv, Q, Vf, Recursion{1e-9, 1e-5})
	fmt.Println(mat.Formatted(Vf))
}
