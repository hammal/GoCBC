package ssm

import (
	"errors"

	"gonum.org/v1/gonum/mat"
)

type AutonomousLinearStateSpaceModel struct {
	// State Dynamics
	A *mat.Dense
	// Observation matrix
	C *mat.Dense
}

// StateObservation returns the observed stateDerivate
// y(t) = C x(t)
// where
// state = x(t) and t is an arbitrary time.
func (model AutonomousLinearStateSpaceModel) StateObservation(t float64, state *mat.VecDense) *mat.VecDense {
	m2, _ := model.A.Dims()
	if m1, _ := state.Dims(); m1 != m2 {
		panic(errors.New("State vector doesn't match state transition matrix"))
	}
	mC, _ := model.C.Dims()
	res := mat.NewVecDense(mC, nil)
	res.MulVec(model.C, state)
	return res
}

// StateDerivative returns the state derivative.
// x'(t) = Ax(t)
// where state = x(t) at an arbitrary time t.
func (model AutonomousLinearStateSpaceModel) StateDerivative(t float64, state *mat.VecDense) *mat.VecDense {
	// Check if state and model parameters match.
	m2, _ := model.A.Dims()
	if m1, _ := state.Dims(); m1 != m2 {
		panic(errors.New("State vector doesn't match state transistion matrix"))
	}

	// Compute state transition
	//  A x(t)
	state.MulVec(model.A, state)

	return state
}
