package control

import (
	"errors"
	"fmt"

	"github.com/hammal/adc/gonumExtensions"
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
	bits []uint
	// Initial state
	state mat.Vector
	// State space model
	StateSpaceModel ssm.StateSpaceModel
	// precomputed control decision vectors for simulation
	controlSimulateLookUp ControlVector
	// precomputed control decision vectors for filtering
	controlFilterLookUpForward  ControlVector
	controlFilterLookUpBackward ControlVector
}

// Simulate the simulation tool for integratorControl
func (c *AnalogSwitchControl) Simulate() [][]float64 {
	var (
		tmpCtrl   mat.Vector
		tmpState  mat.Dense
		tmpSimRes mat.Matrix
	)

	fmt.Println("Starting simulation...")

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
		// fmt.Printf("Control Contribution\n%v\n", mat.Formatted(tmpCtrl))
		// tmpState.Add(tmpState, tmpVec)

		// fmt.Printf("State After \n%v\n", mat.Formatted(&tmpState))

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
	bits := make([]uint, c.NumberOfControls)
	// fmt.Printf("Decisions for \n%v\n => ", mat.Formatted(state))
	for i := 0; i < c.NumberOfControls; i++ {
		tmp = state.AtVec(i) > 0
		// fmt.Printf("%v, ", tmp)
		if tmp {
			bits[i] = 1
		} else {
			bits[i] = 0
		}
	}
	c.bits[index] = bitToIndex(bits)
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

	// tmp := mat.NewVecDense(c.StateSpaceModel.StateSpaceOrder(), nil)

	// retrieve the precomputed vector for 	controlSimulateLookUp [][]mat.Vector [ #control][ decision]
	// var controlIndex, decision int
	// for controlIndex = 0; controlIndex < c.NumberOfControls; controlIndex++ {
	// 	decision = c.bits[index][controlIndex]
	// 	tmp.AddVec(tmp, c.controlSimulateLookUp[controlIndex][decision])
	// }

	tmp := c.controlSimulateLookUp.GetVector(c.bits[index])
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

	tmp := c.controlFilterLookUpForward.GetVector(c.bits[index])

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

	tmp := c.controlFilterLookUpBackward.GetVector(c.bits[index])
	return tmp, nil
}

// GetLength returns the length of control (number of time samples)
func (c AnalogSwitchControl) GetLength() int {
	return len(c.bits)
}

// GetTs returns the sample period
func (c AnalogSwitchControl) GetTs() float64 { return c.Ts }

func (c *AnalogSwitchControl) PreComputeFilterContributions(forwardDynamics, backwardDynamics mat.Matrix) {
	numberOfControlScenarios := (1 << uint(c.NumberOfControls))

	analogswitchForward := analogSwitch{
		systemDynamics: forwardDynamics,
		controls:       c.controls,
		Ts:             c.GetTs(),
	}

	// This is kind of a hack since the differential equation needing solving is
	// dx(t)/dt = -(A + Vb C C^T)x(t) - Bs(t)
	// And the function s(t) was nightmerish since new functions don't get new
	// memory addresses.
	negatedControls := make([]signal.VectorFunction, c.NumberOfControls)
	for index := range negatedControls {
		var tmpVec mat.VecDense
		tmpVec.ScaleVec(-1, c.controls[index].B)
		negatedControls[index] = signal.VectorFunction{
			B: &tmpVec,
			U: c.controls[index].U,
		}
	}

	analogswitchBackward := analogSwitch{
		systemDynamics: backwardDynamics,
		controls:       negatedControls,
		Ts:             c.GetTs(),
	}

	tmpVecForward := make([]mat.Vector, numberOfControlScenarios)
	tmpBoolForward := make([]bool, numberOfControlScenarios)

	tmpVecBackward := make([]mat.Vector, numberOfControlScenarios)
	tmpBoolBackward := make([]bool, numberOfControlScenarios)

	lazycacheForward := lazyCache{
		aSwitch:  analogswitchForward,
		computed: tmpBoolForward,
		cache:    tmpVecForward,
	}

	lazycacheBackward := lazyCache{
		aSwitch:  analogswitchBackward,
		computed: tmpBoolBackward,
		cache:    tmpVecBackward,
	}

	c.controlFilterLookUpForward = &lazycacheForward
	c.controlFilterLookUpBackward = &lazycacheBackward

}

// Returns an initialized analog switch control
func NewAnalogSwitchControl(length int, controls []mat.Vector, ts, t0 float64, state mat.Vector, StateSpaceModel *ssm.LinearStateSpaceModel) *AnalogSwitchControl {
	// IDEA Change controls to reflect resistance values
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

	numberOfPossibleControlCombinations := (2 << uint(numberOfControls))

	// Create decision table
	bits := make([]uint, length)

	analogswitch := analogSwitch{
		systemDynamics: StateSpaceModel.A,
		controls:       ctrl,
		Ts:             ts,
	}

	tmpVec := make([]mat.Vector, numberOfPossibleControlCombinations)
	tmpBool := make([]bool, numberOfPossibleControlCombinations)

	lazycache := lazyCache{
		aSwitch:  analogswitch,
		computed: tmpBool,
		cache:    tmpVec,
	}

	return &AnalogSwitchControl{
		NumberOfControls:      numberOfControls,
		controls:              ctrl,
		Ts:                    ts,
		T0:                    t0,
		bits:                  bits,
		state:                 st,
		StateSpaceModel:       StateSpaceModel,
		controlSimulateLookUp: &lazycache,
	}

}

type analogSwitch struct {
	systemDynamics mat.Matrix
	controls       []signal.VectorFunction
	Ts             float64
}

func (as analogSwitch) GetVector(controlCode uint) mat.Vector {

	ctrlBits := indexToBits(controlCode, len(as.controls))
	ctrlFunction := make([]signal.VectorFunction, len(as.controls))

	// fmt.Printf("Control bits: %v\n", ctrlBits)

	for controlIndex, controlFunction := range as.controls {
		ctrlDecision := (2.*float64(ctrlBits[controlIndex]) - 1.)

		// This is a ugly solution since I can't redefine controlFunction.U
		// Because if so they overlap in memory instead I ended up scaling the whole
		// B vector with the control decision.
		var tmpVec mat.VecDense
		tmpVec.ScaleVec(ctrlDecision, controlFunction.B)

		ctrlFunction[controlIndex] = signal.VectorFunction{
			B: &tmpVec,
			U: controlFunction.U,
		}
	}

	M, _ := as.systemDynamics.Dims()
	linearSystemModel := ssm.NewLinearStateSpaceModel(as.systemDynamics, gonumExtensions.Eye(M, M, 0), ctrlFunction)
	return Solve(linearSystemModel, 0, as.Ts, nil)
}
