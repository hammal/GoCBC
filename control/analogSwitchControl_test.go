package control

import (
	"fmt"
	"testing"

	"github.com/hammal/adc/signal"
	"github.com/hammal/adc/ssm"
	"gonum.org/v1/gonum/mat"
)

func TestAnalogSwitchControl(t *testing.T) {
	order := 4
	length := 10000
	ts := 0.25
	t0 := 0.

	// Create controls
	controls := make([]*mat.VecDense, 4)
	for index, _ := range controls {
		tmp := mat.NewVecDense(order, nil)
		tmp.SetVec(index, 1)
		controls[index] = tmp
	}

	// Create state space Model
	data := make([]float64, order)
	data[0] = 1.
	inp := signal.NewInput(func(arg1 float64) float64 { return 1 }, mat.NewVecDense(order, data))
	stateSpaceModel := ssm.NewIntegratorChain(order, 10, inp)

	control := NewAnalogSwitchControl(length, controls, ts, t0, nil, *stateSpaceModel)

	for index := 0; index < control.NumberOfControls; index++ {
		fmt.Println(mat.Formatted(&control.controlSimulateLookUp[index][0]))
		fmt.Println(mat.Formatted(&control.controlSimulateLookUp[index][1]))
	}

	control.Simulate()

	for index := 0; index < control.Length; index++ {
		for ctrl := 0; ctrl < control.NumberOfControls; ctrl++ {
			fmt.Print(control.bits[index][ctrl])
		}
		fmt.Printf(" for index %v \n", index)
	}

}
