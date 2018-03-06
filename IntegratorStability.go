package main

import (
	"fmt"
	"math"
	"os"
	"strconv"

	"github.com/hammal/adc/control"
	"github.com/hammal/adc/samplingnetwork"
	"gonum.org/v1/gonum/mat"
)

func main() {
	gain, err1 := strconv.ParseFloat(os.Args[1], 64)
	frequency, err2 := strconv.ParseFloat(os.Args[2], 64)
	if err1 != nil && err2 != nil {
		panic("Input not a float")
	}
	fmt.Printf("Simulation For Gain, Frequency;  (%v, %v)\n", gain, frequency)
	phase := math.Pi * 1. / 7.
	resonanceFrequency := frequency
	// oscillator := samplingnetwork.SeriesBlock([]samplingnetwork.SamplingNetwork{samplingnetwork.IntegratorBlock(gain), samplingnetwork.IntegratorBlock(gain)})
	oscillator := samplingnetwork.FeedbackBlock(samplingnetwork.IntegratorBlock(frequency), samplingnetwork.IntegratorBlock(frequency))

	tmpBlock := samplingnetwork.SeriesBlock([]samplingnetwork.SamplingNetwork{samplingnetwork.IntegratorBlock(gain), samplingnetwork.IntegratorBlock(gain)})

	oscillator.System.B = tmpBlock.System.B
	oscillator.Control = tmpBlock.Control

	fmt.Println(mat.Formatted(oscillator.System.A), mat.Formatted(oscillator.System.B), mat.Formatted(oscillator.System.C))
	fmt.Println(oscillator.Control)

	input := []func(float64) float64{
		func(arg float64) float64 {
			return 1. * math.Sin(2*math.Pi*resonanceFrequency*arg+phase)
		},
		// func(arg float64) float64 { return 0. * rand.Float64() },
	}

	length := 100
	controls := oscillator.Control
	ctrl := make([]mat.Vector, len(controls))
	for index := range controls {
		controls[index].SetState(-1.)
		ctrl[index] = controls[index].GetVector()
		fmt.Println(mat.Formatted(ctrl[index]))
	}

	// fmt.Printf("Sanity check, %v, %v, %v, \n", ctrl[0].U(1), ctrl[1].U(1), ctrl[1].U(1))

	ts := 1e-3
	t0 := 0.
	StateSpaceModel := samplingnetwork.LinearSystemToLinearStateSpaceModel(oscillator.System, input)
	// StateSpaceModel.Input[0].B = mat.NewVecDense(2, []float64{0., 0.6 * gain})
	// state := mat.NewVecDense(2, nil)
	oscillatorCtrl := control.NewAnalogSwitchControl(length, ctrl, ts, t0, nil, StateSpaceModel)
	oscillatorCtrl.Simulate()
}
