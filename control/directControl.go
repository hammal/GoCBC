package control

import (
	"errors"

	"github.com/hammal/adc/control/simulator"
	"github.com/hammal/ssm"
	"gonum.org/v1/gonum/mat"
)

type integratorControl struct {
	// Number of samples
	Lenght int
	// Number of controls
	NumberOfControls int
	// Sampling period
	Ts float64
	// Starting time
	T0 float64
	// Bool array with control decisions
	bits [][]bool
	// Simulator
	sim *simulator.Simulator
	// Control signals
	csm []ssm.LinearStateSpaceModel
}

// Simulate the simulation tool for integratorControl
func (c *integratorControl) Simulate() {
	state := mat.NewVecDense(c.csm[0].Order, nil)
	t0 := c.T0
	t1 := t0 + c.Ts
	for index := 0; index < c.Lenght; index++ {
		(*c.sim).Simulate(t0, t1, state, c.csm)
		c.updateControl(state, index)
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
func (c *integratorControl) updateControl(state *mat.VecDense, index int) {
	// Set control bits
	for i := 0; i < c.NumberOfControls; i++ {
		c.bits[index][i] = state.AtVec(i) > 0
	}
	// Update csm
	for i, inp := range c.csm {
		if c.bits[index][i] {
			inp.Input.U = directControlSignalPositive
		} else {
			inp.Input.U = directControlSignalNegative
		}
	}
}

// Copy the bools from internal array and return
func (c *integratorControl) Get(index int) []bool {
	// Check that index exsists
	if index < 0 || index > c.Lenght-1 {
		panic(errors.New("Index out of range"))
	}

	// Allocate new bool array
	tmp := make([]bool, c.NumberOfControls)
	// Populate array by copying from bits
	for bit := 0; bit < c.NumberOfControls; bit++ {
		tmp[bit] = c.bits[index][bit]
	}
	return tmp
}
