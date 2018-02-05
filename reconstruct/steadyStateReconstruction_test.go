package reconstruct

import (
	"fmt"
	"math"
	"testing"

	"github.com/hammal/adc/control"
	"github.com/hammal/adc/signal"
	"github.com/hammal/adc/ssm"
	"gonum.org/v1/gonum/mat"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/plotutil"
	"gonum.org/v1/plot/vg"
)

func TestNewSteadyStateReconstructor(t *testing.T) {
	// runtime.GOMAXPROCS(3)
	N := 7
	beta := 6250.
	b := mat.NewVecDense(N, nil)
	b.SetVec(0, beta)
	input := make([]signal.VectorFunction, 1)
	sig := func(arg1 float64) float64 { return math.Sin(math.Pi*2.*arg1*10. + 0.345) }
	input[0] = signal.NewInput(sig, b)
	sm := ssm.NewIntegratorChain(N, beta, input)

	Length := 500 * 10
	ts := 1. / 16000.
	t0 := 0.
	controls := make([]mat.Vector, N)
	for index := range controls {
		tmp := mat.NewVecDense(N, nil)
		tmp.SetVec(index, -math.Abs(beta))
		controls[index] = tmp
	}

	state := mat.NewVecDense(N, nil)
	ctrl := control.NewAnalogSwitchControl(Length, controls, ts, t0, state, sm)
	fmt.Println("Simulating")
	ctrl.Simulate()

	var inputNoiseCovariance, measurementNoiseCovariance mat.Dense

	sigma_u2 := 1e-8
	sigma_z2 := 1e0
	inputNoiseCovariance.Outer(sigma_u2, b, b)
	measurementNoiseCovariance.Mul(sm.C, sm.C.T())
	measurementNoiseCovariance.Scale(sigma_z2, &measurementNoiseCovariance)

	fmt.Printf("Measurment Noise matrix \n%v\n", mat.Formatted(&measurementNoiseCovariance))
	fmt.Println("PreComputing Reconstruction")
	rec := NewSteadyStateReconstructor(ctrl, &measurementNoiseCovariance, &inputNoiseCovariance, *sm)
	fmt.Println("Reconstructing")
	res := rec.Reconstruction()

	fmt.Println("Forward filter contributions")
	for index := 0; index < ctrl.GetLength(); index++ {
		controlDesc, err := ctrl.GetForwardControlFilterContribution(index)
		// controlDesc, err = ctrl.GetBackwardControlFilterContribution(index)
		if err != nil {
			t.Error(err)
		}
		// for controlIndex, bits := range ctrl.bits[index] {
		// 	fmt.Printf("%v for control %v\n", bits, controlIndex)
		// }
		fmt.Printf("Which results in the contribution\n%v \n", mat.Formatted(controlDesc))

		// fmt.Printf(" for index %v \n", index)
	}

	fmt.Println("Results are")
	for index, term := range res {
		for index2 := range res[index] {
			res[index][index2] *= sigma_u2
		}
		fmt.Printf("S=%v, E=%v\n", sig(float64(index)*ts), term[0])
	}

	sigsig := make([][]float64, Length)
	for index := range sigsig {
		sigsig[index] = []float64{sig(float64(index) * ts)}
	}

	// Plotting
	p, err := plot.New()
	if err != nil {
		panic(err)
	}

	p.Title.Text = "Estimation Comparison for Visual aid"
	p.X.Label.Text = "index"
	p.Y.Label.Text = "."

	err = plotutil.AddLines(p,
		"Estimate", plottify(res)[0],
		"Signal", plottify(sigsig)[0],
	)
	if err != nil {
		panic(err)
	}

	// Save the plot to a PNG file.
	if err := p.Save(4*vg.Inch, 4*vg.Inch, "points.eps"); err != nil {
		panic(err)
	}

}

// randomPoints returns some random x, y points.
func plottify(data [][]float64) []plotter.XYs {
	NumberOfSamples := len(data)
	NumberOfPlots := len(data[0])

	res := make([]plotter.XYs, NumberOfPlots)

	for index := range res {
		pts := make(plotter.XYs, NumberOfSamples)
		for i := range pts {
			pts[i].X = float64(i)
			pts[i].Y = data[i][index]
		}
		res[index] = pts
	}
	return res
}
