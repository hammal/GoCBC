package ssm

import (
	"github.com/hammal/adc/gonumExtensions"
	"github.com/hammal/adc/signal"
	"gonum.org/v1/gonum/mat"
)

type BiLinearStateSpaceModel struct {
	// Linear State Dynamics
	AL mat.Matrix
	// Bi-linear State Dynamics
	AB mat.Matrix
	// Observation matrix
	C mat.Matrix
	// List of input functions
	Input []signal.VectorFunction
}

// Derivative returns the state derivative.
// x'(t) = ALx(t) + Bu(t) + AB Vec(x(t)u(t)^T)
// where state = x(t) at an arbitrary time t. Furthermore, Bu is the input vector field.
func (model BiLinearStateSpaceModel) Derivative(t float64, state mat.Vector) mat.Vector {
	// Define variables
	var (
		tmpState    mat.VecDense
		tmpBilinear mat.VecDense
	)
	// Check if state and model parameters match.
	m2, _ := model.AL.Dims()
	if m1, _ := state.Dims(); m1 != m2 {
		panic("State vector doesn't match state transition matrix")
	}

	// Create Bilinear matrix
	inputVector := mat.NewVecDense(len(model.Input), nil)
	for index := range model.Input {
		inputVector.SetVec(index, model.Input[index].U(t))
	}
	bilinearMatrix := mat.NewDense(m2, len(model.Input), nil)
	// fmt.Printf("State = \n%v\nInputVector = \n%v\n", mat.Formatted(state), mat.Formatted(inputVector))
	bilinearMatrix.Outer(1., state, inputVector)
	bilinearVector := gonumExtensions.Vectorize(bilinearMatrix)
	// fmt.Printf("AB = \n%v\nbilinearVector = \n%v\n", mat.Formatted(model.AB), mat.Formatted(bilinearVector))
	tmpBilinear.MulVec(model.AB, bilinearVector)
	// fmt.Printf("AB = \n%v\nbilinearVector = \n%v\n", mat.Formatted(model.AB), mat.Formatted(bilinearVector))

	// Build input vector
	// tempInput = B input[0](t) ... + B input[N](t)
	tmpInput := mat.NewVecDense(m2, nil)
	for _, input := range model.Input {
		tmpInput.AddVec(tmpInput, input.Bu(t))
	}
	// fmt.Printf("Input contribution \n%v\n", mat.Formatted(tmpInput))

	// Compute state transition
	//  A x(t)
	tmpState.MulVec(model.AL, state)

	// Add to new state derivative vector and return
	tmpInput.AddVec(&tmpState, tmpInput)
	tmpInput.AddVec(&tmpBilinear, tmpInput)

	// fmt.Printf("State contribution \n%v\n", mat.Formatted(&tmpState))
	// fmt.Printf("BiLinear contribution \n%v\n", mat.Formatted(&tmpBilinear))
	// fmt.Println(mat.Formatted(tmpInput))
	return tmpInput
}

// Observation returns the observed stateDerivate
// y(t) = C x(t)
// where
// state = x(t) and t is an arbitrary time.
func (model BiLinearStateSpaceModel) Observation(t float64, state mat.Vector) mat.Vector {
	m2, _ := model.AL.Dims()
	if m1, _ := state.Dims(); m1 != m2 {
		panic("State vector doesn't match state transition matrix")
	}
	mC, _ := model.C.Dims()
	res := mat.NewVecDense(mC, nil)
	res.MulVec(model.C, state)
	return res
}

func (ssm BiLinearStateSpaceModel) StateSpaceOrder() int {
	m, _ := ssm.AL.Dims()
	return m
}

func (ssm BiLinearStateSpaceModel) ObservationSpaceOrder() int {
	m, _ := ssm.C.Dims()
	return m
}

func (ssm BiLinearStateSpaceModel) InputSpaceOrder() int {
	return len(ssm.Input)
}

func (ssm BiLinearStateSpaceModel) Order() int {
	return ssm.StateSpaceOrder()
}
