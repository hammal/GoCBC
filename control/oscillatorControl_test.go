package control

import (
	"math"
	"testing"

	"github.com/hammal/adc/samplingnetwork"
	"github.com/hammal/adc/signal"
)

func TestOscillatorStability(t *testing.T) {
	gain := 1e4
	phase := 0 * math.Pi / 4.
	resonanceFrequency := 2e5
	oscillator := samplingnetwork.OscillatorBlock(gain, resonanceFrequency)

	input := []func(float64) float64{
		func(arg float64) float64 { return 0. / math.Sqrt(2) * math.Sin(2*math.Pi*resonanceFrequency*arg+phase) },
		func(arg float64) float64 { return 0. },
	}

	length := 100
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
	oscillatorCtrl := NewAnalogOscillatorControl(length, ctrl, ts, t0, nil, StateSpaceModel)
	oscillatorCtrl.Simulate()

	// for index := 0; index < oscillatorCtrl.GetLength(); index++ {
	// 	vec, _ := oscillatorCtrl.GetForwardControlFilterContribution(index)
	// 	fmt.Printf("%v\n", mat.Formatted(vec))
	// }
}
