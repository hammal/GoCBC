package control

import (
	"errors"

	"github.com/hammal/adc/ode"
	"github.com/hammal/adc/signal"
	"github.com/hammal/adc/ssm"
	"gonum.org/v1/gonum/mat"
)

// AnalogSwitchControl is the implementation of control with
// open analog switches.
type AnalogSwitchControl struct {
	// Number of controls
	NumberOfControls int
	// Controls
	controls []signal.VectorFunction
	// Sampling period
	Ts float64
	// Starting time
	T0 float64
	// Control decisions
	bits [][]int
	// Initial state
	state mat.Vector
	// State space model
	StateSpaceModel ssm.StateSpaceModel
	// precomputed control decision vectors for simulation
	controlSimulateLookUp [][]mat.Vector
	// precomputed control decision vectors for filtering
	controlFilterLookUpForward  [][]mat.Vector
	controlFilterLookUpBackward [][]mat.Vector
}

// Simulate the simulation tool for integratorControl
func (c *AnalogSwitchControl) Simulate() [][]float64 {
	var (
		tmpCtrl   mat.Vector
		tmpState  mat.Dense
		tmpSimRes mat.Matrix
	)

	res := make([][]float64, c.GetLength())

	tmpState = *mat.NewDense(c.StateSpaceModel.StateSpaceOrder(), 1, nil)

	for row := 0; row < c.StateSpaceModel.StateSpaceOrder(); row++ {
		tmpState.Set(row, 0, c.state.AtVec(row))
	}

	t0 := c.T0
	t1 := t0 + c.Ts
	rk := ode.NewRK4()
	for index := 0; index < c.GetLength(); index++ {
		// fmt.Printf("State Before \n%v\n", mat.Formatted(tmpState))
		// fmt.Printf("Current state = \n%v\n", mat.Formatted(tmpState))
		// Update control based on current state
		c.updateControl(tmpState.ColView(0), index)
		// Simulate the ADC without control
		// There is a complication here since we now the exact state dynamics
		// if the state space model was a linear model. Thus this could be realized
		// using a pre-computed Ad=e^(A Ts) and then using the Runge-Kutta method
		// with zero initial state.
		tmpSimRes, _ = rk.Compute(t0, t1, &tmpState, c.StateSpaceModel)
		// Get the control contributions
		tmpCtrl, _ = c.getControlSimulationContribution(index)
		// Add the control contributions
		// fmt.Printf("Simulation Contribution \n%v\n", mat.Formatted(tmpSimRes))

		// tmpState = mat.NewDense(c.StateSpaceModel.StateSpaceOrder(), 1, nil)
		tmpState.Add(tmpCtrl, tmpSimRes)
		// fmt.Printf("Control Contribution\n%v\n", mat.Formatted(tmpVec))
		// tmpState.Add(tmpState, tmpVec)

		// fmt.Printf("State After \n%v\n", mat.Formatted(tmpState))

		// Move increment to new time step
		t0 += c.Ts
		t1 += c.Ts

		res[index] = make([]float64, c.StateSpaceModel.StateSpaceOrder())
		for row := 0; row < c.StateSpaceModel.StateSpaceOrder(); row++ {
			res[index][row] = tmpState.At(row, 0)
		}
	}
	c.state = tmpState.ColView(0)
	return res
}

// updateControl computes the control decisions for index based on the current
// state, held in the reviver type.
func (c *AnalogSwitchControl) updateControl(state mat.Vector, index int) {
	// Set control bits
	var tmp bool
	// fmt.Printf("Decisions for \n%v\n => ", mat.Formatted(state))
	for i := 0; i < c.NumberOfControls; i++ {
		tmp = state.AtVec(i) > 0
		// fmt.Printf("%v, ", tmp)
		if tmp {
			c.bits[index][i] = 1
		} else {
			c.bits[index][i] = 0
		}

	}
}

// GetControlSimulationContribution returns the control decision vector
// for simulation.
func (c *AnalogSwitchControl) getControlSimulationContribution(index int) (mat.Vector, error) {
	// Check that index exists
	if index < 0 || index > c.GetLength()-1 {
		return nil, errors.New("Index out of range")
	}

	if c.controlSimulateLookUp == nil {
		return nil, errors.New("No pre-computed filter decisions.")
	}

	// fmt.Printf("Number of controls = %v", c.NumberOfControls)

	tmp := mat.NewVecDense(c.StateSpaceModel.StateSpaceOrder(), nil)

	// retrieve the precomputed vector for 	controlSimulateLookUp [][]mat.Vector [ #control][ decision]
	var controlIndex, decision int
	for controlIndex = 0; controlIndex < c.NumberOfControls; controlIndex++ {
		decision = c.bits[index][controlIndex]
		tmp.AddVec(tmp, c.controlSimulateLookUp[controlIndex][decision])
	}
	return tmp, nil
}

func (c AnalogSwitchControl) GetForwardControlFilterContribution(index int) (mat.Vector, error) {
	// Check that index exists
	if index < 0 || index > c.GetLength()-1 {
		return nil, errors.New("index out of range")
	}
	// Check that there are precomputed filter decisions
	if c.controlFilterLookUpForward == nil {
		return nil, errors.New("No pre-computed filter decisions.")
	}

	tmp := mat.NewVecDense(c.StateSpaceModel.StateSpaceOrder(), nil)

	// retrieve the precomputed vector for contolLookUp[ #control][ decision]
	var controlIndex, decision int
	for controlIndex = 0; controlIndex < c.NumberOfControls; controlIndex++ {
		decision = c.bits[index][controlIndex]
		tmp.AddVec(tmp, c.controlFilterLookUpForward[controlIndex][decision])
	}
	return tmp, nil
}

func (c AnalogSwitchControl) GetBackwardControlFilterContribution(index int) (mat.Vector, error) {
	// Check that index exists
	if index < 0 || index > c.GetLength()-1 {
		return nil, errors.New("index out of range")
	}
	// Check that there are precomputed filter decisions
	if c.controlFilterLookUpBackward == nil {
		return nil, errors.New("No pre-computed filter decisions.")
	}

	tmp := mat.NewVecDense(c.StateSpaceModel.StateSpaceOrder(), nil)

	// retrieve the precomputed vector for contolLookUp[ #control][ decision]
	var controlIndex, decision int
	for controlIndex = 0; controlIndex < c.NumberOfControls; controlIndex++ {
		decision = c.bits[index][controlIndex]
		tmp.AddVec(tmp, c.controlFilterLookUpBackward[controlIndex][decision])
	}
	return tmp, nil
}

// GetLength returns the length of control (number of time samples)
func (c AnalogSwitchControl) GetLength() int {
	return len(c.bits)
}

// GetTs returns the sample period
func (c AnalogSwitchControl) GetTs() float64 { return c.Ts }

func (c *AnalogSwitchControl) PreComputeFilterContributions(forwardDynamics, backwardDynamics mat.Matrix) {
	c.controlFilterLookUpForward = make([][]mat.Vector, c.NumberOfControls)
	c.controlFilterLookUpBackward = make([][]mat.Vector, c.NumberOfControls)

	input := make([]signal.VectorFunction, 1)

	var (
		systemForward  *ssm.LinearStateSpaceModel
		systemBackward *ssm.LinearStateSpaceModel
		resFoward      mat.Vector
		resBackward    mat.Vector
	)

	for controlIndex, controlFunction := range c.controls {
		c.controlFilterLookUpForward[controlIndex] = make([]mat.Vector, 2)
		c.controlFilterLookUpBackward[controlIndex] = make([]mat.Vector, 2)
		input[0] = controlFunction

		systemForward = ssm.NewLinearStateSpaceModel(forwardDynamics, forwardDynamics, input)
		resFoward = c.PreCompute(systemForward, 0, c.GetTs())
		negresForward := mat.NewVecDense(c.StateSpaceModel.StateSpaceOrder(), nil)
		negresForward.ScaleVec(-1., resFoward)
		c.controlFilterLookUpForward[controlIndex][0] = negresForward
		c.controlFilterLookUpForward[controlIndex][1] = resFoward

		// fmt.Printf("Forward :\n%v\nNegation: \n%v\n", resFoward, negresForward)

		systemBackward = ssm.NewLinearStateSpaceModel(backwardDynamics, backwardDynamics, input)
		resBackward = c.PreCompute(systemBackward, 0, c.GetTs())
		negresBackward := mat.NewVecDense(c.StateSpaceModel.StateSpaceOrder(), nil)
		negresBackward.ScaleVec(-1, resBackward)
		c.controlFilterLookUpBackward[controlIndex][0] = negresBackward
		c.controlFilterLookUpBackward[controlIndex][1] = resBackward
	}
	// for index, _ := range controlSimulateLookUp {
	// 	// Initialize a receive vector for each control bit
	// 	controlSimulateLookUp[index] = make([]mat.Vector, 2)
}

func (c AnalogSwitchControl) PreCompute(system ode.DifferentiableSystem, from, to float64) mat.Vector {
	o := ode.NewFehlberg45()
	value := mat.NewDense(c.StateSpaceModel.StateSpaceOrder(), 1, nil)
	res, _ := o.Compute(from, to, value, system)
	res2 := res.(*mat.Dense)
	return res2.ColView(0)
}

// Returns an initialized analog switch control
func NewAnalogSwitchControl(length int, controls []mat.Vector, ts, t0 float64, state mat.Vector, StateSpaceModel *ssm.LinearStateSpaceModel) *AnalogSwitchControl {
	order := StateSpaceModel.StateSpaceOrder()
	numberOfControls := len(controls)
	ctrl := make([]signal.VectorFunction, numberOfControls)

	// Construct default controls
	for index := range controls {
		ctrl[index] = signal.NewInput(func(arg1 float64) float64 { return 1. }, controls[index])
	}

	// If state is an nil pointer initialize a new zero vector.
	st, ok := state.(*mat.VecDense)
	if !ok {
		st = mat.NewVecDense(order, nil)
		// state = st
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
	// fmt.Println("Ad: ")
	// fmt.Println(mat.Formatted(Ad))

	// Create decision Simulation lookup table
	controlSimulateLookUp := make([][]mat.Vector, numberOfControls)
	for index, _ := range controlSimulateLookUp {
		// Initialize a receive vector for each control bit
		controlSimulateLookUp[index] = make([]mat.Vector, 2)

		tmp0 := mat.NewVecDense(order, nil)
		tmp1 := mat.NewVecDense(order, nil)

		tmp0.MulVec(Ad, controls[index])
		tmp1.ScaleVec(-1, controls[index])
		tmp1.MulVec(Ad, tmp0)
		controlSimulateLookUp[index][0] = tmp1
		controlSimulateLookUp[index][1] = tmp0

	}

	return &AnalogSwitchControl{
		NumberOfControls:      numberOfControls,
		controls:              ctrl,
		Ts:                    ts,
		T0:                    t0,
		bits:                  bits,
		state:                 st,
		StateSpaceModel:       StateSpaceModel,
		controlSimulateLookUp: controlSimulateLookUp,
	}
}
