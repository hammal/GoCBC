package adc

import (
	"github.com/hammal/adc/control"
	"github.com/hammal/adc/reconstructor"
	"github.com/hammal/adc/simulator"
)

// ADC is the overall interface for controlling the simulation and reconstruction
// of and ADC.
type ADC interface {
	// Simulate the system. Returns a tuple [time index][output]float64
	Simulate() [][]float64
	// Reconstruct the inputs. Returns a tuple [time index][input]float64
	Reconstruct() [][]float64

	// Save the adc system to a file filename
	Save(filename string)
	// Load the adc from a file filename
	Load(filename string) *ADC

	// Get time stamps. Returns time stamps corresponding to the data points from
	// Simulate() and Reconstruct().
	GetTimeStamps() []float64
}

type adc struct {
	cont control.Control
	sim  simulator.Simulator
	rec  reconstructor.Reconstruction
	sys  system
}

// Simulate initates the control for simulating each control sequence
func (a *adc) Simulate() [][]float64 {
	// run the simulation
	a.cont.Simulate()
	// Return the observations
	tmp := a.sim.GetObservations()
	return tmp
}

// Reconstruct using the reconstruction object.
func (a *adc) Reconstruct() [][]float64 {
	return a.rec.Reconstruct()
}

// Return the time stamps for Ts
func (a adc) GetTimeStamps() []float64 {
	N := int((a.sys.t1 - a.sys.t0) / a.sys.ts)
	tmp := make([]float64, N)
	for index := 0; index < N; index++ {
		tmp[index] = float64(index) * a.sys.ts
	}
	return tmp
}
