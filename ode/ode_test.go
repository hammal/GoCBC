package ode

import (
	"fmt"
	"math/rand"
	"testing"

	"github.com/hammal/adc/gonumExtensions"
	"github.com/hammal/adc/signal"
	"github.com/hammal/adc/ssm"
	"gonum.org/v1/gonum/mat"
)

func TestRk4(t *testing.T) {
	test := NewRK4()
	if test.Description.stages != 4 {
		t.Errorf("Not four stages. Rk4 should have four stages. Instead has %v", test.Description.stages)
	}
}

func TestEuler(t *testing.T) {
	test := NewEulerMethod()
	if test.Description.stages != 1 {
		t.Error("Wrong number of stages.")
	}
}

func TestFehlberg(t *testing.T) {
	test := NewFehlberg45()
	if test.Description.stages != 6 {
		t.Error("Wrong number of stages.")
	}
}

func TestCompute(y *testing.T) {
	odeObject := NewRK4()
	t0, t1 := 0., 1.234567891
	initState := mat.NewDense(9, 1, []float64{1.2, 32., 34, 12, 532, 12, 35, 1, 0.91283})
	data := make([]float64, 9*9)
	for index := 0; index < 9*9; index++ {
		data[index] = rand.Float64()
	}
	stateTransitionMatrix := mat.NewDense(9, 9, data)
	inputvector := mat.NewVecDense(9, []float64{1, 0, 0, 0, 0, 0, 0, 0, 0})
	inputs := make([]signal.VectorFunction, 1)
	inputs[0] = signal.NewInput(func(arg1 float64) float64 { return arg1 }, inputvector)
	observationVector := mat.NewDense(1, 9, []float64{0, 0, 0, 0, 0, 0, 0, 0, 1})
	sys := ssm.NewLinearStateSpaceModel(stateTransitionMatrix, observationVector, inputs)
	odeObject.Compute(t0, t1, initState, sys)
	fmt.Println(mat.Formatted(initState))
}

func TestComputeA100States(y *testing.T) {
	odeObject := NewRK4()
	t0, t1 := 0., 1.234567891
	N := 100
	M := 1000
	K := 89
	data := make([]float64, N*M)
	for index := range data {
		data[index] = rand.Float64()
	}
	initState := mat.NewDense(N, M, data)
	data = make([]float64, N*N)
	for index := range data {
		data[index] = rand.Float64()
	}
	inputs := make([]signal.VectorFunction, K)
	for index := range inputs {
		tmpB := mat.NewVecDense(N, nil)
		for index2 := 0; index2 < N; index2++ {
			tmpB.SetVec(index2, rand.ExpFloat64())
		}
		inputs[index] = signal.NewInput(func(arg1 float64) float64 { return arg1 * rand.Float64() }, tmpB)
	}
	sys := ssm.NewIntegratorChain(N, 100.2359, inputs)

	observationVector := gonumExtensions.Eye(M, N, 0)
	sys.C = observationVector
	odeObject.Compute(t0, t1, initState, sys)
	// fmt.Println(mat.Formatted(initState))
}

func TestAdaptiveCompute(y *testing.T) {
	odeObject := NewFehlberg45()
	t0, t1 := 0., 1.
	err := 1e-4
	initState := mat.NewDense(9, 1, []float64{1.2, 32., 34, 12, 532, 12, 35, 1, 0.91283})
	data := make([]float64, 9*9)
	for index := 0; index < 9*9; index++ {
		data[index] = rand.Float64()
	}
	stateTransitionMatrix := mat.NewDense(9, 9, data)
	inputvector := mat.NewVecDense(9, []float64{1, 0, 0, 0, 0, 0, 0, 0, 0})
	inputs := make([]signal.VectorFunction, 1)
	inputs[0] = signal.NewInput(func(arg1 float64) float64 { return 1.24 }, inputvector)

	observationVector := mat.NewDense(1, 9, []float64{0, 0, 0, 0, 0, 0, 0, 0, 1})

	sys := ssm.NewLinearStateSpaceModel(stateTransitionMatrix, observationVector, inputs)

	sucess := odeObject.AdaptiveCompute(t0, t1, err, initState, sys)

	if sucess != nil {
		y.Error(sucess)
	}

	fmt.Println(mat.Formatted(initState))
}

func TestAdaptiveComputeExceedCount(y *testing.T) {
	odeObject := NewFehlberg45()
	t0, t1 := 0., 10.
	err := 1e-9
	initState := mat.NewDense(9, 1, []float64{1.2, 32., 34, 12, 532, 12, 35, 1, 0.91283})
	data := make([]float64, 9*9)
	for index := 0; index < 9*9; index++ {
		data[index] = rand.Float64()
	}
	stateTransitionMatrix := mat.NewDense(9, 9, data)
	inputvector := mat.NewVecDense(9, []float64{1, 0, 0, 0, 0, 0, 0, 0, 0})
	inputs := make([]signal.VectorFunction, 1)
	inputs[0] = signal.NewInput(func(arg1 float64) float64 { return 1.24 }, inputvector)

	observationVector := mat.NewDense(1, 9, []float64{0, 0, 0, 0, 0, 0, 0, 0, 1})

	sys := ssm.NewLinearStateSpaceModel(stateTransitionMatrix, observationVector, inputs)

	sucess := odeObject.AdaptiveCompute(t0, t1, err, initState, sys)

	if sucess == nil {
		y.Error("Didn't produce error for extreme integration length")
	}

	fmt.Println(mat.Formatted(initState))
}

func TestAdaptiveComputeManyStates(y *testing.T) {
	odeObject := NewRK4()
	t0, t1 := 0., 1.234567891
	N := 100
	M := 1000
	K := 89
	err := 1e-9
	data := make([]float64, N*M)
	for index := range data {
		data[index] = rand.Float64()
	}
	initState := mat.NewDense(N, M, data)
	data = make([]float64, N*N)
	for index := range data {
		data[index] = rand.Float64()
	}
	inputs := make([]signal.VectorFunction, K)
	for index := range inputs {
		tmpB := mat.NewVecDense(N, nil)
		for index2 := 0; index2 < N; index2++ {
			tmpB.SetVec(index2, rand.ExpFloat64())
		}
		inputs[index] = signal.NewInput(func(arg1 float64) float64 { return arg1 * rand.Float64() }, tmpB)
	}
	sys := ssm.NewIntegratorChain(N, 100.2359, inputs)

	observationVector := gonumExtensions.Eye(M, N, 0)
	sys.C = observationVector
	odeObject.AdaptiveCompute(t0, t1, err, initState, sys)
	fmt.Println("TestAdaptiveComputeManyStates: Passed")
	// fmt.Println(mat.Formatted(initState))
}
