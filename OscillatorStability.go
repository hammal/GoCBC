package main

import (
	"fmt"
	"math"
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
	gain := 1e2 * (1. + argument)
	fmt.Printf("Simulation For Gain = %v\n", gain)
	phase := math.Pi * 2. / 4.
	resonanceFrequency := 2e5
	oscillator := samplingnetwork.OscillatorBlock(gain, resonanceFrequency)

	input := []func(float64) float64{
		func(arg float64) float64 { return 1 * math.Sin(2*math.Pi*resonanceFrequency*arg+phase) },
		func(arg float64) float64 { return 0. },
	}

	length := 1000
	controls := oscillator.Control
	ctrl := make([]signal.VectorFunction, len(controls))
	for index := range controls {
		controls[index].SetState(-1.)
		ctrl[index] = controls[index].GetResponse()
	}

	ts := 1e-3
	t0 := 0.
	StateSpaceModel := samplingnetwork.LinearSystemToLinearStateSpaceModel(oscillator.System, input)
	// state := mat.NewVecDense(2, nil)
	oscillatorCtrl := control.NewAnalogOscillatorControl(length, ctrl, ts, t0, nil, StateSpaceModel)
	oscillatorCtrl.Simulate()
}
