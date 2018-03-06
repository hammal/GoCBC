package main

import (
	"fmt"
	"math"
	"math/rand"
	"os"
	"strconv"

	"github.com/hammal/adc/control"
	"github.com/hammal/adc/samplingnetwork"
	"github.com/hammal/adc/signal"
)

func main() {
	argument, err := strconv.ParseFloat(os.Args[1], 64)
	if err != nil {
		panic("Input not a float")
	}
	gain := argument
	fmt.Printf("Simulation For Gain = %v\n", gain)
	phase := math.Pi * 1. / 7.
	resonanceFrequency := 2.1 * 1e3
	oscillator := samplingnetwork.OscillatorBlock(gain, resonanceFrequency)

	input := []func(float64) float64{
		func(arg float64) float64 {
			return 8. / argument * math.Sin(2*math.Pi*resonanceFrequency*arg+phase)
		},
		func(arg float64) float64 { return 0.1 / argument * rand.Float64() },
	}

	length := 100000
	controls := oscillator.Control
	ctrl := make([]signal.VectorFunction, len(controls))
	for index := range controls {
		controls[index].SetState(-1.)
		ctrl[index] = controls[index].GetResponse()
	}

	// fmt.Printf("Sanity check, %v, %v, %v, \n", ctrl[0].U(1), ctrl[1].U(1), ctrl[1].U(1))

	ts := 2e-6
	t0 := 0.
	StateSpaceModel := samplingnetwork.LinearSystemToLinearStateSpaceModel(oscillator.System, input)
	// StateSpaceModel.Input[0].B = mat.NewVecDense(2, []float64{0., 0.6 * gain})
	// state := mat.NewVecDense(2, nil)
	oscillatorCtrl := control.NewAnalogOscillatorControl(length, ctrl, ts, t0, nil, StateSpaceModel)
	oscillatorCtrl.Simulate()
}
