package ode

import (
	"fmt"
	"math/rand"
	"testing"

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

func TestCompute(y *testing.T) {
	odeObject := NewRK4()
	t0, t1 := 0., 1.234567891
	initState := mat.NewVecDense(9, []float64{1.2, 32., 34, 12, 532, 12, 35, 1, 0.91283})
	data := make([]float64, 9*9)
	for index := 0; index < 9*9; index++ {
		data[index] = rand.Float64()
	}
	stateTransitionMatrix := mat.NewDense(9, 9, data)
	inputvector := mat.NewVecDense(9, []float64{1, 0, 0, 0, 0, 0, 0, 0, 0})
	inputs := signal.NewInput(func(arg1 float64) float64 { return arg1 }, inputvector)
	observationVector := mat.NewDense(1, 9, []float64{0, 0, 0, 0, 0, 0, 0, 0, 1})
	sys := ssm.NewLinearStateSpaceModel(stateTransitionMatrix, observationVector, inputs)
	odeObject.Compute(t0, t1, initState, sys)
	fmt.Println(mat.Formatted(initState))
}
