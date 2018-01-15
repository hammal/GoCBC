package control

import (
	"errors"

	"github.com/hammal/adc/ode"
	"github.com/hammal/adc/ssm"
	"gonum.org/v1/gonum/mat"
)

type analogSwitchControl struct {
	// Number of samples
	Length int
	// Number of controls
	NumberOfControls int
	// Sampling period
	Ts float64
	// Starting time
	T0 float64
	// Control decisions
	bits [][]int
	// Intial state
	state *mat.VecDense
	// State space model
	stateSpaceModel ssm.StateSpaceModel
	// precomputed control descision vectors
	controlLookUp [][]mat.Vector
}

// Simulate the simulation tool for integratorControl
func (c *analogSwitchControl) Simulate() {
	t0 := c.T0
	t1 := t0 + c.Ts
	rk := ode.NewRK4()
	var tmpCtrl []mat.Vector
	for index := 0; index < c.Length; index++ {
		// Update control based on current state
		c.updateControl(index)
		// Simulate the ADC without control
		rk.Compute(t0, t1, c.state, c.stateSpaceModel)
		// Get the control contributions
		c.GetControlContribution(index)
		// Add the control contributions
		for _, tmpVec := range tmpCtrl {
			c.state.AddVec(c.state, tmpVec)
		}
		// Move increment to new time step
		t0 += c.Ts
		t1 += c.Ts
	}
}

// These are the two direct control zero order hold
// filters.
func directControlSignalPositive(t float64) float64 {
	return 1.
}

func directControlSignalNegative(t float64) float64 {
	return -1.
}

// This is the control logic for the integratorControl
func (c *analogSwitchControl) updateControl(index int) {
	// Set control bits
	var tmp bool
	for i := 0; i < c.NumberOfControls; i++ {
		tmp = c.state.AtVec(i) > 0
		if tmp {
			c.bits[index][i] = 1
		} else {
			c.bits[index][i] = 0
		}

	}
}

// Copy the bools from internal array and return
func (c *analogSwitchControl) GetControlContribution(index int) []*mat.Vector {
	// Check that index exsists
	if index < 0 || index > c.Length-1 {
		panic(errors.New("Index out of range"))
	}

	tmp := make([]*mat.Vector, c.NumberOfControls)

	// retrive the precomputed vector for contolLookUp[ #control][ decision]
	var controlIndex, decision int
	for controlIndex = 0; controlIndex < c.NumberOfControls; controlIndex++ {
		decision = c.bits[index][controlIndex]
		tmp[controlIndex] = &c.controlLookUp[controlIndex][decision]
	}
	return tmp
}
