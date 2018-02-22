package control

import (
	"errors"
	"fmt"
	"math"

	"github.com/hammal/adc/gonumExtensions"
	"github.com/hammal/adc/ode"
	"github.com/hammal/adc/signal"
	"github.com/hammal/adc/ssm"
	"gonum.org/v1/gonum/mat"
)

type oscillatorControl struct {
	U func(float64) float64
	B mat.Vector
	C mat.Vector
}

// AnalogSwitchControl is the implementation of control with
// open analog switches.
type OscillatingControl struct {
	// Number of controls
	NumberOfControls int
	// Controls
	controls []oscillatorControl
	// Sampling period
	Ts float64
	// Starting time
	T0 float64
	// Control decisions
	bits []uint
	// Initial state
	State mat.Vector
	// State space model
	StateSpaceModel ssm.BiLinearStateSpaceModel
	// Inputs
	Inputs []signal.VectorFunction
	// FilterContributions
	controlFilterLookUpForward  oscillatorSwitch
	controlFilterLookUpBackward oscillatorSwitch
}

// Simulate the simulation tool for integratorControl
func (c *OscillatingControl) Simulate() [][]float64 {
	var (
		tmpState  mat.Dense
		tmpSimRes mat.Matrix
	)

	fmt.Println("Starting simulation...")

	res := make([][]float64, c.GetLength())

	tmpState = *mat.NewDense(c.StateSpaceModel.StateSpaceOrder(), 1, nil)

	for row := 0; row < c.StateSpaceModel.StateSpaceOrder(); row++ {
		tmpState.Set(row, 0, c.State.AtVec(row))
	}

	t0 := c.T0
	t1 := t0 + c.Ts
	var err error
	// rk := ode.NewRK4()
	rk := ode.NewFehlberg45()
	for index := 0; index < c.GetLength(); index++ {
		// fmt.Printf("State Before \n%v\n", mat.Formatted(tmpState))
		// fmt.Printf("Current state = \n%v\n", mat.Formatted(tmpState))
		// Update control based on current state
		c.updateControl(tmpState.ColView(0), index)
		tmpCtrl := make([]signal.VectorFunction, c.NumberOfControls+c.StateSpaceModel.InputSpaceOrder())
		tmpCtrl2, _ := c.getControlSimulationContribution(index)
		tmpCtrl = append(c.Inputs, tmpCtrl2...)
		c.StateSpaceModel.Input = tmpCtrl

		// fmt.Printf("tmpState Before = \n%v\n", mat.Formatted(&tmpState))
		// tmpSimRes, _ = rk.Compute(t0, t1, &tmpState, c.StateSpaceModel)
		tmpSimRes, err = rk.AdaptiveCompute(t0, t1, 1e-3, &tmpState, c.StateSpaceModel)
		if err != nil {
			panic(err)
		}
		// fmt.Printf("tmpSimRes \n%v\n", mat.Formatted(tmpSimRes))
		M, N := tmpSimRes.Dims()
		for row := 0; row < M; row++ {
			for column := 0; column < N; column++ {
				tmpState.Set(row, column, tmpSimRes.At(row, column))
			}
		}

		t0 += c.Ts
		t1 += c.Ts

		res[index] = make([]float64, c.StateSpaceModel.StateSpaceOrder())
		for row := 0; row < c.StateSpaceModel.StateSpaceOrder(); row++ {
			res[index][row] = tmpSimRes.At(row, 0)
		}
	}

	c.State = tmpState.ColView(0)
	return res
}

// updateControl computes the control decisions for index based on the current
// state, held in the reviver type.
func (c *OscillatingControl) updateControl(state mat.Vector, index int) {
	// Set control bits
	var tmp bool
	var tmpFloat float64
	I := gonumExtensions.Eye(c.StateSpaceModel.Order(), c.StateSpaceModel.Order(), 0)
	bits := make([]uint, c.NumberOfControls)

	// fmt.Println()
	for i := 0; i < c.NumberOfControls; i++ {
		tmpFloat = mat.Inner(c.controls[i].C, I, state)
		tmp = tmpFloat > 0.
		// fmt.Printf("%v, %v ", tmpFloat, tmp)
		if tmp {
			bits[i] = 1
		} else {
			bits[i] = 0
		}
	}
	// fmt.Print("\n")
	// fmt.Println(bits[0], bits[1])
	c.bits[index] = bitToIndex(bits)

	// These are all for monitoring controllability.

	a3 := mat.Inner(c.controls[0].C, I, state)
	a4 := mat.Inner(c.controls[1].C, I, state)
	e1 := math.Pow(state.AtVec(0), 2)
	e2 := math.Pow(state.AtVec(1), 2)
	phase := math.Atan2(state.AtVec(1), state.AtVec(0)) / math.Pi * 180
	phaseError := math.Atan2(a3, a4) / math.Pi * 180.

	fmt.Printf("\nEnergy: Total = %5.e, PerState = (%5.e, %5.e,) [V^2], Amplitude of Buffers = (%+2.3e, %+2.3e) [V] and Phase / PhaseError = %+4.f / %+4.f [deg]", e1+e2, e1, e2, a3, a4, phase, phaseError)

	if index > 0 {
		if c.bits[index] != c.bits[index-1] {
			fmt.Print("\tControl Switch!\t")
			bits1 := indexToBits(c.bits[index-1], len(bits))
			for index2 := range bits {
				if bits1[index2] != bits[index2] {
					fmt.Printf("Nr %d,\t", index2)
				}
			}
		}
	}
}

// getControlSimulationContribution returns the control decision vector
// for simulation.
func (c *OscillatingControl) getControlSimulationContribution(index int) ([]signal.VectorFunction, error) {
	// Check that index exists
	if index < 0 || index > c.GetLength()-1 {
		return nil, errors.New("Index out of range")
	}

	ctrlBits := indexToBits(c.bits[index], len(c.controls))
	ctrlFunction := make([]signal.VectorFunction, len(c.controls))

	for controlIndex, controlFunction := range c.controls {
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

	return ctrlFunction, nil
}

func (c OscillatingControl) GetForwardControlFilterContribution(index int) (mat.Vector, error) {
	// Check that index exists
	if index < 0 || index > c.GetLength()-1 {
		return nil, errors.New("index out of range")
	}

	// tmp := c.controlFilterLookUpForward.GetVector(c.bits[index])
	t0 := c.Ts*(float64(index-1)) + c.T0
	tmp := c.controlFilterLookUpForward.GetVector(c.bits[index], t0, t0+c.Ts)

	return tmp, nil
}

func (c OscillatingControl) GetBackwardControlFilterContribution(index int) (mat.Vector, error) {
	// Check that index exists
	if index < 0 || index > c.GetLength()-1 {
		return nil, errors.New("index out of range")
	}

	t0 := c.Ts*(float64(index)) + c.T0
	tmp := c.controlFilterLookUpForward.GetVector(c.bits[index], t0, t0+c.Ts)
	return tmp, nil
}

// GetLength returns the length of control (number of time samples)
func (c OscillatingControl) GetLength() int {
	return len(c.bits)
}

// GetTs returns the sample period
func (c OscillatingControl) GetTs() float64 { return c.Ts }

func (c *OscillatingControl) PreComputeFilterContributions(forwardDynamics, backwardDynamics mat.Matrix) {

	oscillatorSwitchForward := oscillatorSwitch{
		systemDynamics: forwardDynamics,
		controls:       oscillatorSwitchToAnalogSwitch(c.controls),
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

	ocillatorSwitchBackward := oscillatorSwitch{
		systemDynamics: backwardDynamics,
		controls:       negatedControls,
		Ts:             c.GetTs(),
	}

	c.controlFilterLookUpForward = oscillatorSwitchForward
	c.controlFilterLookUpBackward = ocillatorSwitchBackward

}

// Returns an initialized analog switch control
func NewAnalogOscillatorControl(length int, controls []signal.VectorFunction, ts, t0 float64, state mat.Vector, StateSpaceModel *ssm.LinearStateSpaceModel) *OscillatingControl {
	if state != nil {
		panic("This is not implemented. Pass a nil value instead")
	}

	order := StateSpaceModel.StateSpaceOrder()
	numberOfControls := len(controls)
	ctrl := make([]oscillatorControl, numberOfControls)

	// Construct default controls
	for index := range controls {
		tmpC := mat.NewVecDense(order+numberOfControls, nil)
		tmpC.SetVec(index+order, 1.)

		tmpB := mat.NewVecDense(order+numberOfControls, nil)
		for row := 0; row < order; row++ {
			tmpB.SetVec(row, controls[index].B.AtVec(row))
		}
		ctrl[index] = oscillatorControl{
			U: controls[index].U,
			B: tmpB,
			C: tmpC,
		}
	}

	// New state space models
	var tmp1, tmp2, ALnew mat.Dense
	tmp1.Augment(StateSpaceModel.A, mat.NewDense(order, numberOfControls, nil))
	// ALnew.Stack(&tmp1, mat.NewDense(numberOfControls, numberOfControls+order, nil))

	data := make([]float64, numberOfControls)
	for index := range data {
		data[index] = -1000000.
	}
	// This is to define leaky integrators --> low pass filters for the bilinear terms
	I := mat.NewDiagonal(numberOfControls, data)
	tmp2.Augment(mat.NewDense(numberOfControls, order, nil), I)
	ALnew.Stack(&tmp1, &tmp2)

	ABnew := mat.NewDense(order+numberOfControls, (order+numberOfControls)*(numberOfControls+len(StateSpaceModel.Input)), nil)

	for controlIndex := range controls {
		for candidateState := 0; candidateState < order+numberOfControls; candidateState++ {
			// Find all states where the control is added
			if ctrl[controlIndex].B.AtVec(candidateState) != 0 {
				// TODO there is several issues with this definition. works only for one oscillator block
				ABnew.Set(order+numberOfControls-1-controlIndex, (order+numberOfControls)*(controlIndex+len(StateSpaceModel.Input))+candidateState, 1.)
			}
		}
	}

	fmt.Printf("New oscillator with additional states AL, AB \n%v\n%v\n", mat.Formatted(&ALnew), mat.Formatted(ABnew))

	// Adjust inputs
	tmpInput := make([]signal.VectorFunction, len(StateSpaceModel.Input))
	for index := range StateSpaceModel.Input {
		tmpB := mat.NewVecDense(order+numberOfControls, nil)
		for row := 0; row < order; row++ {
			tmpB.SetVec(row, StateSpaceModel.Input[index].B.AtVec(row))
		}
		tmpInput[index] = signal.VectorFunction{
			U: StateSpaceModel.Input[index].U,
			B: tmpB,
		}
	}

	// Create a new State Space Model

	NewStateSpaceModel := ssm.BiLinearStateSpaceModel{
		AL:    &ALnew,
		AB:    ABnew,
		Input: tmpInput,
		C:     gonumExtensions.Eye(order+numberOfControls, order+numberOfControls, 0),
	}

	// If state is an nil pointer initialize a new zero vector.
	var st *mat.VecDense
	st = mat.NewVecDense(order+numberOfControls, nil)

	// Create decision table
	bits := make([]uint, length)

	return &OscillatingControl{
		NumberOfControls: numberOfControls,
		controls:         ctrl,
		Ts:               ts,
		T0:               t0,
		bits:             bits,
		State:            st,
		StateSpaceModel:  NewStateSpaceModel,
		Inputs:           tmpInput,
	}

}

//
type oscillatorSwitch struct {
	systemDynamics mat.Matrix
	controls       []signal.VectorFunction
	Ts             float64
}

func (as oscillatorSwitch) GetVector(controlCode uint, from, to float64) mat.Vector {

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
	return Solve(linearSystemModel, from, to, nil)
}

func oscillatorSwitchToAnalogSwitch(controls []oscillatorControl) []signal.VectorFunction {
	res := make([]signal.VectorFunction, len(controls))
	for index, control := range controls {
		res[index] = signal.VectorFunction{
			B: control.B,
			U: control.U,
		}
	}
	return res
}
