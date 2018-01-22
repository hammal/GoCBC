package signal

import (
	"gonum.org/v1/gonum/mat"
)

// VectorFunction is an input stateSpaceModel abstraction that instead of returning
// a scalar associates each argument with a vector valued output. For instance
// in the state space model:
// X'(t) = AX(t) + BU(t)
// BU(t) is a vectorial function decomposed as a scalar function U(t)-> Reals
// and a vector B \in Reals^N. Where N is the dimension of the state space.
// DU(t) as with BU(t) is the vectorial function decomposition for the state Observation.
type VectorFunction struct {
	U func(float64) float64
	B mat.Vector
}

// Bu is an alias in the state space model:
//
// Ax + Bu(t)
// where Bu(t) is a vectorial function
func (vf VectorFunction) Bu(t float64) mat.Vector {
	return vf.Value(t)
}

// Value returns the vectorial function value
func (vf VectorFunction) Value(t float64) mat.Vector {
	var res mat.VecDense
	res.CloneVec(vf.B)
	res.ScaleVec(vf.U(t), &res)
	return &res
}

// func (vf VectorFunction) Du(t float64) *mat.VecDense {
// 	m, _ := vf.D.Dims()
// 	res := mat.NewVecDense(m, nil)
// 	res.CloneVec(vf.D)
// 	res.ScaleVec(vf.u(t), res)
// 	return res
// }

// NewInput returns a pointer to a new VectorFunction object
// initalised with u(t) and B
func NewInput(u func(float64) float64, B mat.Vector) VectorFunction {
	return VectorFunction{u, B}
}
