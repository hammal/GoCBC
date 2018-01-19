package ssm

import (
	"errors"
	"sync"

	"github.com/hammal/adc/signal"
	"gonum.org/v1/gonum/mat"
)

// LinearStateSpaceModel struct represent the system
//
// x'(t) = A x(t) + B input[0](t) ... + B input[N](t)
//
// y(t) = C x(t)
//
// where N is the number of inputs
type LinearStateSpaceModel struct {
	// State Dynamics
	A mat.Matrix
	// Observation matrix
	C mat.Matrix
	// List of input functions
	Input []signal.VectorFunction
}

// NewIntegratorChain returns a linear state space model of an integrator chain
// of size N with input.
func NewIntegratorChain(N int, stageGain float64, input []signal.VectorFunction) *LinearStateSpaceModel {
	a := make([]float64, N*N)
	c := make([]float64, N)
	stride := N
	for row := 0; row < N; row++ {
		c[row] = 1
		for column := 0; column < N; column++ {
			if row == (column + 1) {
				a[row*stride+column] = stageGain
			}
		}
	}
	A := mat.NewDense(N, N, a)
	C := mat.NewDense(1, N, c)
	return NewLinearStateSpaceModel(A, C, input)
}

// NewLinearStateSpaceModel creates a new Linear state space model
func NewLinearStateSpaceModel(A, C mat.Matrix, input []signal.VectorFunction) *LinearStateSpaceModel {
	// Check that system parameters match
	m, n := A.Dims()
	_, nC := C.Dims()
	mB, _ := input[0].B.Dims()
	if m != n || nC != m || mB != m {
		panic(errors.New("System Parameters don't match"))
	}
	// Check that input dimensions match

	sys := LinearStateSpaceModel{A, C, input}
	return &sys
}

// Derivative returns the state derivative.
// x'(t) = Ax(t) + Bu(t)
// where state = x(t) at an arbitrary time t. Furthermore, Bu is the input vector field.
func (model LinearStateSpaceModel) Derivative(t float64, state mat.Vector) mat.Vector {
	// Define variables
	var (
		tmpState mat.VecDense
		tmpInput mat.VecDense
	)

	// Check if state and model parameters match.
	m2, _ := model.A.Dims()
	if m1, _ := state.Dims(); m1 != m2 {
		panic(errors.New("State vector doesn't match state transition matrix"))
	}

	// Build input vector
	// tempInput = B input[0](t) ... + B input[N](t)
	tmpInput = *mat.NewVecDense(m2, nil)
	for _, input := range model.Input {
		tmpInput.AddVec(&tmpInput, input.Bu(t))
	}
	// Compute state transition
	//  A x(t)
	tmpState.MulVec(model.A, state)

	// Add to new state derivative vector and return
	tmpInput.AddVec(&tmpState, &tmpInput)
	return &tmpInput
}

// Observation returns the observed stateDerivate
// y(t) = C x(t)
// where
// state = x(t) and t is an arbitrary time.
func (model LinearStateSpaceModel) Observation(t float64, state mat.Vector) mat.Vector {
	m2, _ := model.A.Dims()
	if m1, _ := state.Dims(); m1 != m2 {
		panic(errors.New("State vector doesn't match state transition matrix"))
	}
	mC, _ := model.C.Dims()
	res := mat.NewVecDense(mC, nil)
	res.MulVec(model.C, state)
	return res
}

// computeStateTransistion computes the e^(At) where A is a square matrix and
// t is a scalar.
func computeStateTransistion(t float64, m *mat.Dense) {
	m.Scale(t, m)
	m.Exp(m)
}

// ImpulseResponse computes the impulse response of the system and returns it in
// an array [numberOfObservations][numberOfInputs][tap at time t]float64
func (model LinearStateSpaceModel) ImpulseResponse(t []float64) [][][]float64 {
	var wg sync.WaitGroup
	numberOfTaps := len(t)
	numberOfObservations, _ := model.C.Dims()
	numberOfInputs := 1

	// Initalise an 3D array --> (numberOfObservations X numberOfInputs X numberOfTaps)
	res := make([][][]float64, numberOfObservations)
	for obs := range res {
		res[obs] = make([][]float64, numberOfInputs)
		for inp := range res[obs] {
			res[obs][inp] = make([]float64, numberOfTaps)
		}
	}

	wg.Add(numberOfTaps)
	for index, time := range t {
		// Compute the different impulses as a go routine
		go func(i int, t float64) {
			defer wg.Done()
			var tmp, tmp2 mat.Dense
			tmp.Clone(model.A)
			// tmp = e^(tmp * t)
			computeStateTransistion(t, &tmp)
			// C e^(tmp * t)
			tmp2.Mul(model.C, &tmp)
			var tempInput mat.VecDense
			for _, input := range model.Input {
				tempInput.MulVec(&tmp2, input.B)
			}
			for obs := range res {
				res[obs][0][i] = tempInput.AtVec(obs)
			}
		}(index, time)
	}
	wg.Wait()
	return res
}

func (ssm LinearStateSpaceModel) StateSpaceOrder() int {
	m, _ := ssm.A.Dims()
	return m
}

func (ssm LinearStateSpaceModel) ObservationSpaceOrder() int {
	m, _ := ssm.C.Dims()
	return m
}

func (ssm LinearStateSpaceModel) InputSpaceOrder() int {
	return len(ssm.Input)
}

// func (model LinearStateSpaceModel) FrequencyResponse(f []float64) [][][]complex128 {
// 	numberOfTaps := len(f)
// 	numberOfObservations := model.Order
// 	numberOfInputs := len(model.Input)
//
// 	// Initalise an 3D array --> (numberOfObservations X numberOfInputs X numberOfTaps)
// 	res := make([][][]complex128, numberOfObservations)
// 	for obs := range res {
// 		res[obs] = make([][]complex128, numberOfInputs)
// 		for inp := range res[obs] {
// 			res[obs][inp] = make([]complex128, numberOfTaps)
// 		}
// 	}

// Have not implemented this since I don't fully understand how to work out
// complex matrices
// panic(errors.New("Not yet implemented"))

// for index, frequency := range f {
// 	// Compute the different frequency response as a go routine
// 	go func(i int, f float64) {
// 		tmp := mat.NewDense(model.order, model.order, nil)
// 		tmp.Clone(model.A)
// 		// tmp = e^(tmp * t)
// 		computeStateTransistion(t, tmp)
// 		// C e^(tmp * t)
// 		tmp.Mul(model.C, tmp)
// 		for inp, inputObject := range model.Input {
// 			tempInput := mat.NewVecDense(model.order, nil)
// 			tempInput.MulVec(tmp, inputObject.B)
// 			tempInput.MulVec(model.C, tempInput)
// 			for obs := range res {
// 				res[obs][inp][index] = tempInput.AtVec(obs)
// 			}
// 		}
// 	}(index, frequency)
// }
// return res

// }
