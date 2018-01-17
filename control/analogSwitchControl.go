package control

import (
	"errors"

	"github.com/hammal/adc/ode"
	"github.com/hammal/adc/ssm"
	"gonum.org/v1/gonum/mat"
)

// Returns an initialized analog switch control
func NewAnalogSwitchControl(length int, controls []*mat.VecDense, ts, t0 float64, state *mat.VecDense, StateSpaceModel ssm.LinearStateSpaceModel) *AnalogSwitchControl {
	order := StateSpaceModel.StateSpaceOrder()
	numberOfControls := len(controls)

	// If state is an nil pointer initialize a new zero vector.
	if state == nil {
		state = mat.NewVecDense(order, nil)
	}

	// Create decision table
	bits := make([][]int, length)
	for bitIndex, _ := range bits {
		bits[bitIndex] = make([]int, numberOfControls)
	}

	// Decision Simulation Lookup table
	// var tmp0, tmp1 *mat.VecDense
	// Compute e^(A Ts)
	Ad := zeroOrderHold(StateSpaceModel.A, ts)

	// Create decision Simulation lookup table
	controlSimulateLookUp := make([][]mat.VecDense, numberOfControls)
	for index, _ := range controlSimulateLookUp {
		// Initialize a receive vector for each control bit
		controlSimulateLookUp[index] = make([]mat.VecDense, 2)

		tmp0 := mat.NewVecDense(order, nil)
		tmp1 := mat.NewVecDense(order, nil)

		tmp1.MulVec(Ad, controls[index])
		tmp0.ScaleVec(-1, controls[index])
		tmp0.MulVec(Ad, tmp0)
		controlSimulateLookUp[index][0] = *tmp0
		controlSimulateLookUp[index][1] = *tmp1

	}

	// ControplFilterLookUp
	controlFilterLookUp := make([][]mat.VecDense, numberOfControls)
	for index, _ := range controlFilterLookUp {
		controlFilterLookUp[index] = make([]mat.VecDense, 2)
	}

	return &AnalogSwitchControl{
		numberOfControls,
		ts,
		t0,
		bits,
		state,
		StateSpaceModel,
		controlSimulateLookUp,
		controlFilterLookUp,
	}
}

// AnalogSwitchControl is the implementation of control with
// open analog switches.
type AnalogSwitchControl struct {
	// Number of controls
	NumberOfControls int
	// Sampling period
	Ts float64
	// Starting time
	T0 float64
	// Control decisions
	bits [][]int
	// Initial state
	state *mat.VecDense
	// State space model
	StateSpaceModel ssm.StateSpaceModel
	// precomputed control descision vectors for simulation
	controlSimulateLookUp [][]mat.VecDense
	// precomputed control descision vectors for filtering
	controlFilterLookUp [][]mat.VecDense
}

// Simulate the simulation tool for integratorControl
func (c *AnalogSwitchControl) Simulate() {
	t0 := c.T0
	t1 := t0 + c.Ts
	rk := ode.NewRK4()
	var tmpCtrl []*mat.VecDense
	for index := 0; index < c.GetLength(); index++ {
		// Update control based on current state
		c.updateControl(index)
		// Simulate the ADC without control
		// There is a complication here since we now the exact state dynamics
		// if the state space model was a linear model. Thus this could be realized
		// using a pre-computed Ad=e^(A Ts) and then using the Runge-Kutta method
		// with zero initial state.
		rk.Compute(t0, t1, c.state, c.StateSpaceModel)
		// Get the control contributions
		tmpCtrl = c.getControlSimulationContribution(index)
		// Add the control contributions
		for _, tmpVec := range tmpCtrl {
			c.state.AddVec(c.state, tmpVec)
		}
		// Move increment to new time step
		t0 += c.Ts
		t1 += c.Ts
	}
}

// updateControl computes the control decisions for index based on the current
// state, held in the reviver type.
func (c *AnalogSwitchControl) updateControl(index int) {
	// Set control bits
	var tmp bool
	for i := 0; i < c.NumberOfControls; i++ {
		tmp = c.state.AtVec(i) > 0
		if tmp {
			c.bits[index][i] = 0
		} else {
			c.bits[index][i] = 1
		}

	}
}

// GetControlSimulationContribution returns the control decision vector
// for simulation.
func (c *AnalogSwitchControl) getControlSimulationContribution(index int) []*mat.VecDense {
	// Check that index exists
	if index < 0 || index > c.GetLength()-1 {
		panic(errors.New("Index out of range"))
	}

	tmp := make([]*mat.VecDense, c.NumberOfControls)

	// retrieve the precomputed vector for contolLookUp[ #control][ decision]
	var controlIndex, decision int
	for controlIndex = 0; controlIndex < c.NumberOfControls; controlIndex++ {
		decision = c.bits[index][controlIndex]
		tmp[controlIndex] = &c.controlSimulateLookUp[controlIndex][decision]
	}
	return tmp
}

// GetControlFilterContribution returns the control decision vector
// for simulation.
func (c *AnalogSwitchControl) GetControlFilterContribution(index int) ([]*mat.VecDense, error) {
	// Check that index exists
	if index < 0 || index > c.GetLength()-1 {
		return nil, errors.New("index out of range")
	}
	// Check that there are precomputed filter decisions
	if c.controlFilterLookUp == nil {
		return nil, errors.New("No pre-computed filter decisions.")
	}

	tmp := make([]*mat.VecDense, c.NumberOfControls)

	// retrieve the precomputed vector for contolLookUp[ #control][ decision]
	var controlIndex, decision int
	for controlIndex = 0; controlIndex < c.NumberOfControls; controlIndex++ {
		decision = c.bits[index][controlIndex]
		tmp[controlIndex] = &c.controlFilterLookUp[controlIndex][decision]
	}
	return tmp, nil
}

// GetLength returns the length of control (number of time samples)
func (c AnalogSwitchControl) GetLength() int {
	return len(c.bits)
}

// GetTs returns the sample period
func (c AnalogSwitchControl) GetTs() float64 { return c.Ts }

func (c AnalogSwitchControl) GetForwardControlFilterContribution(index int) []mat.Vector {
	// TODO
	panic("Not yet implemented")
	return nil
}

func (c AnalogSwitchControl) GetBackwardControlFilterContribution(index int) []mat.Vector {
	// TODO
	panic("Not yet implemented")
	return nil
}

func (c *AnalogSwitchControl) PreComputeFilterContributions(forwardDynamics, backwardDynamics mat.Matrix) {
	// TODO
	panic("Not yet implementd")
}
