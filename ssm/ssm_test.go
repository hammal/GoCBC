package ssm

import (
	"errors"
	"fmt"
	"testing"

	"github.com/hammal/adc/signal"
	"gonum.org/v1/gonum/mat"
)

func TestNewIntegratorChain(t *testing.T) {
	N := 5
	data := make([]float64, N)
	data[0] = 1.
	inp := signal.NewInput(func(arg1 float64) float64 { return 1 }, mat.NewVecDense(N, data))
	stateSpaceModel := NewIntegratorChain(N, 10, inp)
	fmt.Print(mat.Formatted(stateSpaceModel.A))
	var zero mat.Dense
	zero.Pow(stateSpaceModel.A, N)
	for row := 0; row < N; row++ {
		for col := 0; col < N; col++ {
			if zero.At(row, col) != 0 {
				fmt.Print(mat.Formatted(&zero))
				panic(errors.New("Not an integrator chain"))
			}
		}
	}
}

func TestStateSpaceModel(t *testing.T) {
	A := mat.NewDense(3, 3, []float64{0, 0, 0, 1, 0, 0, 0, 1, 0})
	B := mat.NewVecDense(3, []float64{1, 0, 0})
	C := mat.NewDense(1, 3, []float64{0, 0, 1})

	fmt.Printf("State Space Model:\nA = \n%v\nB = \n%v\nC = \n%v\n", mat.Formatted(A), mat.Formatted(B), mat.Formatted(C))

	inputs := signal.NewInput(func(t float64) float64 { return float64(t) }, B)

	// Test LinearStateSpaceModel
	stateSpaceModel := NewLinearStateSpaceModel(A, C, inputs)
	state := mat.NewVecDense(3, nil)
	for i := 0; i < 100; i++ {
		sD := stateSpaceModel.StateDerivative(float64(i), state)
		state.AddVec(sD, state)
		sO := stateSpaceModel.StateObservation(float64(i), state)
		// sO := 5
		fmt.Printf("State = %v\n, for time %v\nStateDerivative = %v\nStateObservation = %v\n", mat.Formatted(state, mat.Prefix("        ")), i, mat.Formatted(sD, mat.Prefix("                  ")), mat.Formatted(sO))
	}
}

func TestImpulseResponse(t *testing.T) {
	fmt.Println("\n\n\nTestImpulseResponse")
	A := mat.NewDense(3, 3, []float64{0, 0, 0, 1, 0, 0, 0, 1, 0})
	B := mat.NewVecDense(3, []float64{1, 0, 0})
	C := mat.NewDense(1, 3, []float64{0, 0, 1})

	fmt.Printf("State Space Model:\nA = \n%v\nB = \n%v\nC = \n%v\n", mat.Formatted(A), mat.Formatted(B), mat.Formatted(C))

	inputs := signal.NewInput(func(t float64) float64 { return float64(t) }, B)

	// Test LinearStateSpaceModel
	stateSpaceModel := NewLinearStateSpaceModel(A, C, inputs)
	time := make([]float64, 10000)
	var Ts float64
	Ts = 0.002319379123781927312987
	for index := range time {
		time[index] = float64(index) * Ts
	}
	impulseResponse := stateSpaceModel.ImpulseResponse(time)
	fmt.Print(impulseResponse)
}

func BenchmarkImpulseResponse(b *testing.B) {
	A := mat.NewDense(3, 3, []float64{1, 0, 0, 1, 0, 0, 0, 1, 0})
	B := mat.NewVecDense(3, []float64{1, 0, 0})
	C := mat.NewDense(1, 3, []float64{0, 0, 1})
	inputs := signal.NewInput(func(t float64) float64 { return float64(t) }, B)

	// Test LinearStateSpaceModel
	stateSpaceModel := NewLinearStateSpaceModel(A, C, inputs)
	time := make([]float64, 10000)
	var Ts float64
	Ts = 0.002319379123781927312987
	for index := range time {
		time[index] = float64(index) * Ts
	}
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		stateSpaceModel.ImpulseResponse(time)
	}
}
