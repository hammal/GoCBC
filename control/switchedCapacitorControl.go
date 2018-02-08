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

var R float64 = 100.
var C float64 = 100e-6

// AnalogSwitchControl is the implementation of control with
// open analog switches.
type SwitchedCapacitorControl struct {
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
	StateSpaceModel *ssm.LinearStateSpaceModel
	// precomputed control decision vectors for simulation
	controlSimulateLookUp ControlVector
	// precomputed control decision vectors for filtering
	controlFilterLookUpForward  ControlVector
	controlFilterLookUpBackward ControlVector
}

// Simulate the simulation tool for integratorControl
func (c *SwitchedCapacitorControl) Simulate() [][]float64 {
	var (
		tmpCtrl   mat.Vector
		tmpState  mat.Dense
		tmpSimRes mat.Matrix
	)

	fmt.Println("Starting simulation...")

	res := make([][]float64, c.GetLength())

	fmt.Println(c.StateSpaceModel.StateSpaceOrder())
	fmt.Println(mat.Formatted(c.StateSpaceModel.A))
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

		fmt.Printf("State After \n%v\n", mat.Formatted(&tmpState))

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
func (c *SwitchedCapacitorControl) updateControl(state mat.Vector, index int) {
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
func (c *SwitchedCapacitorControl) getControlSimulationContribution(index int) (mat.Vector, error) {
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

func (c SwitchedCapacitorControl) GetForwardControlFilterContribution(index int) (mat.Vector, error) {
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

func (c SwitchedCapacitorControl) GetBackwardControlFilterContribution(index int) (mat.Vector, error) {
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
func (c SwitchedCapacitorControl) GetLength() int {
	return len(c.bits)
}

// GetTs returns the sample period
func (c SwitchedCapacitorControl) GetTs() float64 { return c.Ts }

func (c *SwitchedCapacitorControl) PreComputeFilterContributions(forwardDynamics, backwardDynamics mat.Matrix) {
	numberOfControlScenarios := (2 << uint(c.NumberOfControls))

	analogswitchForward := capacativeSwitch{
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

	analogswitchBackward := capacativeSwitch{
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
func NewSwitchedCapacitorControl(length int, controls []SwitchedCapacitor, ts, t0 float64, state mat.Vector, StateSpaceModel *ssm.LinearStateSpaceModel) *SwitchedCapacitorControl {
	order := StateSpaceModel.StateSpaceOrder()
	numberOfControls := len(controls)
	ctrl := make([]signal.VectorFunction, numberOfControls)
	// inputs := make([]signal.VectorFunction, len(StateSpaceModel.Input))

	for index := range controls {
		ctrl[index] = signal.VectorFunction{
			U: func(float64) float64 { return 0 },
			B: controls[index].B,
		}
	}

	// If state is an nil pointer initialize a new zero vector.
	st, ok := state.(*mat.VecDense)
	if !ok {
		st = mat.NewVecDense(order+numberOfControls, nil)
		// state = st
	}

	numberOfPossibleControlCombinations := (1 << uint(numberOfControls))

	// Create decision table
	bits := make([]uint, length)

	analogswitch := capacativeSwitch{
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

	return &SwitchedCapacitorControl{
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

// SwitchedCapacitor describes the SC circuit
// associated with each control. The R and C value
// Sets the decay rate i.e.
// dV/dt = - V/RC
// Furthermore, the B steering vector additionally considers the
// 1/R C_integrator of the integrator.
type SwitchedCapacitor struct {
	R, C float64
	B    mat.Vector
}

type capacativeSwitch struct {
	systemDynamics mat.Matrix
	controls       []signal.VectorFunction
	Ts             float64
}

func (as capacativeSwitch) GetVector(controlCode uint) mat.Vector {

	ctrlBits := indexToBits(controlCode, len(as.controls))

	order, _ := as.systemDynamics.Dims()
	numberOfControls := len(as.controls)

	ctrl := make([]signal.VectorFunction, numberOfControls)

	// Construct additional states for the SC memory
	scValues := make([]float64, numberOfControls)
	scControlConnection := mat.NewDense(order, numberOfControls, nil)

	for index := range scValues {
		scValues[index] = -1. / (C * R)
		for row := 0; row < order; row++ {
			scControlConnection.Set(row, index, as.controls[index].B.AtVec(row))
		}
	}
	scStates := mat.NewDiagonal(numberOfControls, scValues)

	var Anew, tmpMat1, tmpMat2 mat.Dense
	tmpMat1.Augment(as.systemDynamics, scControlConnection)
	tmpMat2.Augment(mat.NewDense(numberOfControls, order, nil), scStates)
	Anew.Stack(&tmpMat1, &tmpMat2)
	// var Cnew mat.Dense
	// tmpMat1.Reset()
	// tmpMat1.Augment(StateSpaceModel.C, mat.NewDense(order, numberOfControls, nil))
	// tmpMat2.Scale(0, &tmpMat2)
	// Cnew.Stack(&tmpMat1, &tmpMat2)
	// Cnew := gonumExtensions.Eye(order+numberOfControls, order+numberOfControls, 0)

	// fmt.Println(mat.Formatted(&Anew))
	// fmt.Println(mat.Formatted(&Cnew))

	// Construct default controls
	for index := range as.controls {
		tmpVec := mat.NewVecDense(order+numberOfControls, nil)
		for row := 0; row < numberOfControls; row++ {
			tmpVec.SetVec(row+order, as.controls[index].B.AtVec(row))
		}
		ctrl[index] = signal.NewInput(func(arg1 float64) float64 { return 0. }, tmpVec)
	}

	// Adjust inputs
	// for index := range StateSpaceModel.Input {
	// 	tmpVec := mat.NewVecDense(order+numberOfControls, nil)
	// 	for row := 0; row < order; row++ {
	// 		tmpVec.SetVec(row, as.controls[index].B.AtVec(row))
	// 	}
	// 	inputs[index] = signal.VectorFunction{
	// 		U: as.controls[index].U,
	// 		B: tmpVec,
	// 	}
	// }

	// NewSSM := ssm.NewLinearStateSpaceModel(&Anew, Cnew, inputs)

	controlState := mat.NewDense(order+numberOfControls, 1, nil)

	for controlIndex, controlFunction := range ctrl {
		ctrlDecision := (2.*float64(ctrlBits[controlIndex]) - 1.)

		// This is a ugly solution since I can't redefine controlFunction.U
		// Because if so they overlap in memory instead I ended up scaling the whole
		// B vector with the control decision.
		var tmpVec mat.VecDense
		tmpVec.ScaleVec(ctrlDecision, unityVector(controlFunction.B))

		controlState.Add(controlState, &tmpVec)
	}
	dummyInput := make([]signal.VectorFunction, 1)
	dummyInput[0] = signal.VectorFunction{
		B: ctrl[0].B,
		U: func(arg float64) float64 { return 0. },
	}

	linearSystemModel := ssm.NewLinearStateSpaceModel(&Anew, gonumExtensions.Eye(order+numberOfControls, order+numberOfControls, 0), dummyInput)

	tmpRes := Solve(linearSystemModel, 0, as.Ts, controlState)
	// Now as this is an approximate solution the states assosiated with then
	// capacitors needs to be drained as they would if implemented as an SC-circuit.
	res := mat.NewVecDense(order+numberOfControls, nil)
	for index := 0; index < order; index++ {
		res.SetVec(index, tmpRes.AtVec(index))
	}
	return res.SliceVec(0, order)
}

func unityVector(vector mat.Vector) mat.Vector {
	var tmpVec mat.VecDense
	tmpVec.MulElemVec(vector, vector)
	norm := mat.Norm(&tmpVec, 1)
	tmpVec.ScaleVec(1./norm, &tmpVec)
	return &tmpVec
}
