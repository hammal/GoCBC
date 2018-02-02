package control

import (
	"fmt"
	"math"
	"testing"

	"github.com/hammal/adc/signal"
	"github.com/hammal/adc/ssm"
	"gonum.org/v1/gonum/mat"
)

func TestAnalogSwitchControlSimulate(t *testing.T) {
	order := 4
	length := 1000
	ts := 1. / 16000.
	t0 := 0.

	// Create controls
	controls := make([]mat.Vector, 4)
	for index, _ := range controls {
		tmp := mat.NewVecDense(order, nil)
		tmp.SetVec(index, -6250.)
		controls[index] = tmp
	}

	// Create state space Model
	data := make([]float64, order)
	data[0] = -6250.
	inp := make([]signal.VectorFunction, 1)
	inp[0] = signal.NewInput(func(arg1 float64) float64 { return math.Exp(-arg1) }, mat.NewVecDense(order, data))
	// inp[0] = signal.NewInput(func(arg1 float64) float64 { return 0. }, mat.NewVecDense(order, data))
	stateSpaceModel := ssm.NewIntegratorChain(order, -6250, inp)

	ctrl := NewAnalogSwitchControl(length, controls, ts, t0, nil, stateSpaceModel)

	for index := 0; index < ctrl.NumberOfControls; index++ {
		fmt.Println(mat.Formatted(ctrl.controlSimulateLookUp[index][0]))
		fmt.Println(mat.Formatted(ctrl.controlSimulateLookUp[index][1]))
	}

	ctrl.Simulate()

	// for index := 0; index < control.Length(); index++ {
	// 	for ctrl := 0; ctrl < control.NumberOfControls; ctrl++ {
	// 		fmt.Print(control.bits[index][ctrl])
	// 	}
	// 	fmt.Printf(" for index %v \n", index)
	// }

}

func TestAnalogSwitchControlPreComputeFilter(t *testing.T) {
	order := 4
	length := 100
	ts := 1. / 16000.
	t0 := 0.

	// Create controls
	controls := make([]mat.Vector, 4)
	for index, _ := range controls {
		tmp := mat.NewVecDense(order, nil)
		tmp.SetVec(index, -6250.)
		controls[index] = tmp
	}

	// Create state space Model
	data := make([]float64, order)
	data[0] = -6250.
	inp := make([]signal.VectorFunction, 1)
	inp[0] = signal.NewInput(func(arg1 float64) float64 { return math.Exp(-arg1) }, mat.NewVecDense(order, data))
	// inp[0] = signal.NewInput(func(arg1 float64) float64 { return 0. }, mat.NewVecDense(order, data))
	stateSpaceModel := ssm.NewIntegratorChain(order, -6250, inp)

	ctrl := NewAnalogSwitchControl(length, controls, ts, t0, nil, stateSpaceModel)

	ctrl.PreComputeFilterContributions(stateSpaceModel.A, stateSpaceModel.A)

	for index := 0; index < ctrl.NumberOfControls; index++ {
		fmt.Println("Forward: ")
		fmt.Println(mat.Formatted(ctrl.controlFilterLookUpForward[index][0]))
		fmt.Println(mat.Formatted(ctrl.controlFilterLookUpForward[index][1]))
		fmt.Println("Backward: ")
		fmt.Println(mat.Formatted(ctrl.controlFilterLookUpBackward[index][0]))
		fmt.Println(mat.Formatted(ctrl.controlFilterLookUpBackward[index][1]))
	}

}

func TestAnalogSwitchControlSimulateAndGetFilterContributions(t *testing.T) {
	order := 4
	length := 5
	ts := 1. / 16000.
	t0 := 0.

	// Create controls
	controls := make([]mat.Vector, 4)
	for index, _ := range controls {
		tmp := mat.NewVecDense(order, nil)
		tmp.SetVec(index, -6250.)
		controls[index] = tmp
	}

	// Create state space Model
	data := make([]float64, order)
	data[0] = -6250.
	inp := make([]signal.VectorFunction, 1)
	inp[0] = signal.NewInput(func(arg1 float64) float64 { return math.Exp(-arg1) }, mat.NewVecDense(order, data))
	// inp[0] = signal.NewInput(func(arg1 float64) float64 { return 0. }, mat.NewVecDense(order, data))
	stateSpaceModel := ssm.NewIntegratorChain(order, -6250, inp)

	ctrl := NewAnalogSwitchControl(length, controls, ts, t0, nil, stateSpaceModel)
	//
	// for index := 0; index < ctrl.NumberOfControls; index++ {
	// 	fmt.Println(mat.Formatted(ctrl.controlSimulateLookUp[index][0]))
	// 	fmt.Println(mat.Formatted(ctrl.controlSimulateLookUp[index][1]))
	// }

	ctrl.Simulate()
	ctrl.PreComputeFilterContributions(stateSpaceModel.A, stateSpaceModel.A)

	fmt.Println("Forward filter contributions")
	for index := 0; index < ctrl.GetLength(); index++ {
		controlDesc, err := ctrl.GetForwardControlFilterContribution(index)
		controlDesc, err = ctrl.GetBackwardControlFilterContribution(index)
		if err != nil {
			t.Error(err)
		}
		for controlIndex, bits := range ctrl.bits[index] {
			fmt.Printf("%v for control %v\n", bits, controlIndex)
		}
		fmt.Printf("Which results in the contribution\n%v \n", mat.Formatted(controlDesc))

		// fmt.Printf(" for index %v \n", index)
	}
}
