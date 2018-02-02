package reconstruct

import (
	"fmt"
	"math"
	"testing"

	"github.com/hammal/adc/control"
	"github.com/hammal/adc/signal"
	"github.com/hammal/adc/ssm"
	"gonum.org/v1/gonum/mat"
)

func TestNewSteadyStateReconstructor(t *testing.T) {
	// runtime.GOMAXPROCS(3)
	N := 5
	beta := -6250.
	b := mat.NewVecDense(N, nil)
	b.SetVec(0, beta)
	input := make([]signal.VectorFunction, 1)
	input[0] = signal.NewInput(func(arg1 float64) float64 { return math.Sin(math.Pi*2.*arg1 + 0.345) }, b)
	sm := ssm.NewIntegratorChain(N, beta, input)

	Length := 10000000
	ts := 1. / 16000.
	t0 := 0.
	controls := make([]mat.Vector, N)
	for index, _ := range controls {
		tmp := mat.NewVecDense(N, nil)
		tmp.SetVec(index, -math.Abs(beta))
		controls[index] = tmp
	}

	state := mat.NewVecDense(N, nil)
	ctrl := control.NewAnalogSwitchControl(Length, controls, ts, t0, state, sm)
	fmt.Println("Simulating")
	ctrl.Simulate()

	var inputNoiseCovariance, measurementNoiseCovariance mat.Dense

	sigma_u2 := 1e-3
	sigma_z2 := 1e-2
	inputNoiseCovariance.Outer(sigma_u2, b, b)
	measurementNoiseCovariance.Mul(sm.C, sm.C.T())
	measurementNoiseCovariance.Scale(sigma_z2, &measurementNoiseCovariance)

	fmt.Printf("Measurment Noise matrix \n%v\n", mat.Formatted(&measurementNoiseCovariance))
	fmt.Println("PreComputing Reconstruction")
	rec := NewSteadyStateReconstructor(ctrl, &measurementNoiseCovariance, &inputNoiseCovariance, *sm)
	fmt.Println("Reconstructing")
	rec.Reconstruction()
}
