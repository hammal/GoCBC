package ssm

import (
	"math"

	"gonum.org/v1/gonum/mat"
)

// VectorFunction is an input ssm abstraction that instead of returning
// a scalar assoisates each argument with a vector valued output. For instance
// in the state space model:
// X'(t) = AX(t) + BU(t)
// BU(t) is a vectorial function decomposed as a scalar function U(t)-> Reals
// and a vector B \in Reals^N. Where N is the dimension of the state space.
// DU(t) as with BU(t) is the vectorial function decomposition for the state Observation.
type VectorFunction struct {
	U func(float64) float64
	B *mat.VecDense
}

// Bu returns the vectorial function
// B u(t)
func (vf VectorFunction) Bu(t float64) *mat.VecDense {
	m, _ := vf.B.Dims()
	res := mat.NewVecDense(m, nil)
	res.CloneVec(vf.B)
	res.ScaleVec(vf.U(t), res)
	return res
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
func NewInput(u func(float64) float64, B *mat.VecDense) VectorFunction {
	return VectorFunction{u, B}
}

// DiracDelta is a dirac delta distrubution as defined in
// https://en.wikipedia.org/wiki/Dirac_delta_function
func DiracDelta(x float64) float64 {
	// These could all be done offline
	var a = 1e-9
	var a2 = a * a
	var C1 = 1. / (math.Abs(a) * math.Sqrt(math.Pi))
	// Return distribution value at x
	return C1 * math.Exp(-x*x/a2)
}
