package control

import (
	"fmt"
	"math"
	"testing"

	"github.com/hammal/adc/signal"
	"github.com/hammal/adc/ssm"
	"gonum.org/v1/gonum/mat"
)

func TestSwitchedCapacitorSimulate(t *testing.T) {
	order := 4
	length := 1000
	ts := 1. / 16000.
	t0 := 0.

	// Create controls
	controls := make([]SwitchedCapacitor, 4)
	for index, _ := range controls {
		tmp := mat.NewVecDense(order, nil)
		tmp.SetVec(index, -6250.)
		controls[index] = SwitchedCapacitor{
			R: 100.,
			C: 100e-6,
			B: tmp,
		}
	}

	// Create state space Model
	data := make([]float64, order)
	data[0] = -6250.
	inp := make([]signal.VectorFunction, 1)
	inp[0] = signal.NewInput(func(arg1 float64) float64 { return math.Exp(-arg1) }, mat.NewVecDense(order, data))
	// inp[0] = signal.NewInput(func(arg1 float64) float64 { return 0. }, mat.NewVecDense(order, data))
	stateSpaceModel := ssm.NewIntegratorChain(order, -6250, inp)

	ctrl := NewSwitchedCapacitorControl(length, controls, ts, t0, nil, stateSpaceModel)

	for index := 0; index < (1 << uint(ctrl.NumberOfControls)); index++ {
		fmt.Printf("Control Response to index = %v, bitcode %v\n", index, indexToBits(uint(index), 4))
		fmt.Println(mat.Formatted(ctrl.controlSimulateLookUp.GetVector(uint(index))))
	}

	ctrl.Simulate()

	for index := 0; index < ctrl.GetLength(); index++ {
		fmt.Printf("Control codeword %v for index %v\n", ctrl.bits[index], index)
	}

}

// func TestAnalogSwitchControlPreComputeFilter(t *testing.T) {
// 	order := 4
// 	length := 100
// 	ts := 1. / 16000.
// 	t0 := 0.
//
// 	// Create controls
// 	controls := make([]mat.Vector, 4)
// 	for index, _ := range controls {
// 		tmp := mat.NewVecDense(order, nil)
// 		tmp.SetVec(index, -6250.)
// 		controls[index] = tmp
// 	}
//
// 	// Create state space Model
// 	data := make([]float64, order)
// 	data[0] = -6250.
// 	inp := make([]signal.VectorFunction, 1)
// 	inp[0] = signal.NewInput(func(arg1 float64) float64 { return math.Exp(-arg1) }, mat.NewVecDense(order, data))
// 	// inp[0] = signal.NewInput(func(arg1 float64) float64 { return 0. }, mat.NewVecDense(order, data))
// 	stateSpaceModel := ssm.NewIntegratorChain(order, -6250, inp)
//
// 	ctrl := NewAnalogSwitchControl(length, controls, ts, t0, nil, stateSpaceModel)
//
// 	ctrl.PreComputeFilterContributions(stateSpaceModel.A, stateSpaceModel.A)
//
// 	for index := 0; index < (2 << uint(ctrl.NumberOfControls)); index++ {
// 		fmt.Println("Forward: ")
// 		fmt.Println(mat.Formatted(ctrl.controlFilterLookUpForward.GetVector(uint(index))))
// 		fmt.Println("Backward: ")
// 		fmt.Println(mat.Formatted(ctrl.controlFilterLookUpBackward.GetVector(uint(index))))
// 	}
//
// }
//
// func TestAnalogSwitchControlSimulateAndGetFilterContributions(t *testing.T) {
// 	order := 4
// 	length := 5
// 	ts := 1. / 16000.
// 	t0 := 0.
//
// 	// Create controls
// 	controls := make([]mat.Vector, 4)
// 	for index, _ := range controls {
// 		tmp := mat.NewVecDense(order, nil)
// 		tmp.SetVec(index, -6250.)
// 		controls[index] = tmp
// 	}
//
// 	// Create state space Model
// 	data := make([]float64, order)
// 	data[0] = -6250.
// 	inp := make([]signal.VectorFunction, 1)
// 	inp[0] = signal.NewInput(func(arg1 float64) float64 { return math.Exp(-arg1) }, mat.NewVecDense(order, data))
// 	// inp[0] = signal.NewInput(func(arg1 float64) float64 { return 0. }, mat.NewVecDense(order, data))
// 	stateSpaceModel := ssm.NewIntegratorChain(order, -6250, inp)
//
// 	ctrl := NewAnalogSwitchControl(length, controls, ts, t0, nil, stateSpaceModel)
// 	//
// 	// for index := 0; index < ctrl.NumberOfControls; index++ {
// 	// 	fmt.Println(mat.Formatted(ctrl.controlSimulateLookUp[index][0]))
// 	// 	fmt.Println(mat.Formatted(ctrl.controlSimulateLookUp[index][1]))
// 	// }
//
// 	ctrl.Simulate()
// 	ctrl.PreComputeFilterContributions(stateSpaceModel.A, stateSpaceModel.A)
//
// 	fmt.Println("Forward filter contributions")
// 	for index := 0; index < ctrl.GetLength(); index++ {
// 		controlDesc, err := ctrl.GetForwardControlFilterContribution(index)
// 		controlDesc, err = ctrl.GetBackwardControlFilterContribution(index)
// 		if err != nil {
// 			t.Error(err)
// 		}
// 		// for controlIndex, bits := range ctrl.bits[index] {
// 		// 	fmt.Printf("%v for control %v\n", bits, controlIndex)
// 		// }
// 		fmt.Printf("Control %v results in the contribution\n%v \n", indexToBits(uint(index), order), mat.Formatted(controlDesc))
//
// 		// fmt.Printf(" for index %v \n", index)
// 	}
// }
//
// func TestBitsToIndex(t *testing.T) {
// 	indices := []uint{1, 2, 3, 4, 5, 6}
// 	bits := [][]uint{{1}, {0, 1}, {1, 1}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1}}
//
// 	for index := range indices {
// 		tmp := bitToIndex(bits[uint(index)])
// 		fmt.Printf("tmp = %v\n", tmp)
// 		if indices[index] != tmp {
// 			t.Error(fmt.Sprintf("%v is not equal to %v\\n", indices[index], tmp))
// 		}
//
// 		tmp2 := indexToBits(indices[index], len(bits[index]))
// 		fmt.Printf("tmp2 = %v\n", tmp2)
// 		for index2 := range tmp2 {
// 			if tmp2[index2] != bits[index][index2] {
// 				t.Error(fmt.Sprintf("%v is not equal to %v\\n", tmp2, bits[index]))
// 			}
// 		}
// 	}
//
// 	// for index := 100; index < 1000; index++ {
// 	// 	if index != bitToIndex(indexToBits(index, int(math.Log2(float64(index))))) {
// 	// 		t.Error("Not reversible")
// 	// 	}
// 	// }
// }
