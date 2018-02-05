package reconstruct

import (
	"fmt"
	"math/rand"
	"testing"
	"time"

	"github.com/hammal/adc/gonumExtensions"
	"github.com/hammal/adc/signal"
	"github.com/hammal/adc/ssm"
	"gonum.org/v1/gonum/mat"
)

// func TestCareRecursive(t *testing.T) {
// 	N := 4
// 	stageGain := 10.
// 	sigma_u := 1e-5
// 	sigma_z := 1e-3
// 	B := mat.NewVecDense(N, []float64{1, 0, 0, 0})
// 	input := make([]signal.VectorFunction, 1)
// 	input[0] = signal.NewInput(func(arg1 float64) float64 {
// 		return 1.
// 	}, B)
// 	stateSpaceModel := ssm.NewIntegratorChain(N, stageGain, input)
// 	Q := mat.NewDense(N, N, nil)
// 	Q.Outer(sigma_u, B, B)
//
// 	var Rinv mat.Dense
// 	Rinv.Scale(1./sigma_z, gonumExtensions.Eye(N, N, 0))
//
// 	Vf := mat.NewDense(N, N, nil)
// 	for index := 0; index < N; index++ {
// 		Vf.Set(index, index, 1.)
// 	}
// 	var Res mat.Matrix
// 	var tmp1 mat.Dense
// 	var R mat.Dense
//
// 	tmp1.Mul(&Rinv, stateSpaceModel.C)
// 	R.Mul(stateSpaceModel.C.T(), &tmp1)
// 	Res = care(stateSpaceModel.A.T(), &R, Q, Vf, Recursion{1e-3, 1e-5})
// 	fmt.Println(mat.Formatted(Res))
// }

func TestCareEigen(t *testing.T) {
	rand.Seed(time.Now().UTC().UnixNano())
	N := 2
	stageGain := 10.
	sigma_u := 1e-5
	sigma_z := 1e-3
	B := mat.NewVecDense(N, []float64{1, 0})
	input := make([]signal.VectorFunction, 1)
	input[0] = signal.NewInput(func(arg1 float64) float64 {
		return 1.
	}, B)
	stateSpaceModel := ssm.NewIntegratorChain(N, stageGain, input)
	data := make([]float64, N*N)
	for i := range data {
		data[i] = rand.NormFloat64()
	}
	// stateSpaceModel.A = mat.NewDense(N, N, data)

	fmt.Println(mat.Formatted(stateSpaceModel.A))
	Q := mat.NewDense(N, N, nil)
	Q.Outer(sigma_u, B, B)

	fmt.Printf("Q = \n%v\n", mat.Formatted(Q))

	var Rinv mat.Dense
	Rinv.Scale(1./sigma_z, gonumExtensions.Eye(N, N, 0))

	Vf := mat.NewDense(N, N, nil)
	for index := 0; index < N; index++ {
		Vf.Set(index, index, 1.)
	}
	var Res mat.Matrix
	var tmp1 mat.Dense
	var R mat.Dense

	tmp1.Mul(&Rinv, stateSpaceModel.C)
	R.Mul(stateSpaceModel.C.T(), &tmp1)

	fmt.Printf("C = \n%v\nR^(-1) = \n%v\nR = \n%v\n", mat.Formatted(stateSpaceModel.C), mat.Formatted(&Rinv), mat.Formatted(&R))

	Res = care(stateSpaceModel.A.T(), &R, Q, Vf, MatrixFactorization{})
	fmt.Println(mat.Formatted(Res))
}

func TestFiniteDuration(t *testing.T) {
	N := 4
	stageGain := 10.
	sigma_u := 1e-5
	sigma_z := 1e-3
	B := mat.NewVecDense(N, []float64{1, 0, 0, 0})
	input := make([]signal.VectorFunction, 1)
	input[0] = signal.NewInput(func(arg1 float64) float64 {
		return 1.
	}, B)
	stateSpaceModel := ssm.NewIntegratorChain(N, stageGain, input)
	Q := mat.NewDense(N, N, nil)
	Q.Outer(sigma_u, B, B)

	var Rinv mat.Dense
	Rinv.Scale(1./sigma_z, gonumExtensions.Eye(N, N, 0))

	Vf := mat.NewDense(N, N, nil)
	for index := 0; index < N; index++ {
		Vf.Set(index, index, 1.)
	}
	var Res mat.Matrix
	var tmp1 mat.Dense
	var R mat.Dense

	tmp1.Mul(&Rinv, stateSpaceModel.C)
	R.Mul(stateSpaceModel.C.T(), &tmp1)

	// Res = care(stateSpaceModel.A.T(), &R, Q, Vf, FiniteDuration{})
	fmt.Println(mat.Formatted(Res))
}

func TestForDavisionAndMakiPaper(t *testing.T) {
	N := 9
	dataA := [81]float64{
		0, 1, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0.2165, -0.0356, 0, -0.0299, 0, -0.027, 0,
		-0.458, 1, -0.0133, 0.0004, 0, 0.0006, 0, 0.0007, 0,
		0, 0, 0, 0, 1, 0, 0, 0, 0,
		0, 0, 0, -29.81, -0.0546, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 1, 0, 0,
		0, 0, 0, 0, 0, -169, -0.13, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 1,
		0, 0, 0, 0, 0, 0, 0, -334.3, -0.1828,
	}

	dataB := [9]float64{
		0, -1.138, -0.0348, 0, 29.56, 0, 47.25, 0, 16.40,
	}

	dataQ := [9]float64{
		0.1, 0.05, 0.5, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4,
	}

	dataR := [1]float64{1.}

	B := mat.NewVecDense(N, dataB[:])
	input := make([]signal.VectorFunction, 1)
	input[0] = signal.NewInput(func(arg1 float64) float64 {
		return 1.
	}, B)
	stateSpaceModel := ssm.NewLinearStateSpaceModel(mat.NewDense(N, N, dataA[:]), gonumExtensions.Eye(N, N, 0), input)

	fmt.Printf("A = \n%v\nC = \n%v\n", mat.Formatted(stateSpaceModel.A), mat.Formatted(stateSpaceModel.C))

	Q := mat.NewDiagonalRect(N, N, dataQ[:])

	fmt.Printf("Q = \n%v\n", mat.Formatted(Q))

	Rinv := mat.NewDense(1, 1, dataR[:])

	Vf := mat.NewDense(N, N, nil)
	for index := 0; index < N; index++ {
		Vf.Set(index, index, 1.)
	}
	var Res mat.Matrix
	var tmp1 mat.Dense
	var R mat.Dense

	tmp1.Mul(Rinv, stateSpaceModel.Input[0].B.T())
	R.Mul(stateSpaceModel.Input[0].B, &tmp1)

	fmt.Printf("R^(-1) = \n%v\nR = \n%v\n", mat.Formatted(Rinv), mat.Formatted(&R))

	Res = care(stateSpaceModel.A.T(), &R, Q, Vf, MatrixFactorization{})
	fmt.Println(mat.Formatted(Res))
}

// func TestCareNewtonMethod(t *testing.T) {
// 	N := 2
// 	stageGain := 1000.
// 	sigma_u := 1e-3
// 	sigma_z := 1e-3
// 	vec := make([]float64, N)
// 	vec[0] = stageGain
// 	B := mat.NewVecDense(N, vec)
// 	input := make([]signal.VectorFunction, 1)
// 	input[0] = signal.NewInput(func(arg1 float64) float64 {
// 		return 1.
// 	}, B)
// 	stateSpaceModel := ssm.NewIntegratorChain(N, stageGain, input)
// 	Q := mat.NewDense(N, N, nil)
// 	Q.Outer(sigma_u, B, B)
//
// 	var Rinv mat.Dense
// 	Rinv.Scale(1./sigma_z, gonumExtensions.Eye(N, N, 0))
//
// 	Vf := mat.NewDense(N, N, nil)
// 	for index := 0; index < N; index++ {
// 		Vf.Set(index, index, 10000)
// 	}
// 	var Res mat.Matrix
// 	var tmp1 mat.Dense
// 	var R mat.Dense
//
// 	tmp1.Mul(&Rinv, stateSpaceModel.C)
// 	R.Mul(stateSpaceModel.C.T(), &tmp1)
// 	// Vf = care(stateSpaceModel.A.T(), &R, Q, Vf, MatrixFactorization{})
// 	Res = care(stateSpaceModel.A.T(), &R, Q, Vf, NewtonMethod{precision: 1e-10})
// 	fmt.Println(mat.Formatted(Res))
// }

func TestCompareMathematica(t *testing.T) {
	N := 2
	stageGain := 1000.
	sigma_u := 1e-6
	sigma_z := 1e-0
	vec := make([]float64, N)
	vec[0] = stageGain
	B := mat.NewVecDense(N, vec)
	input := make([]signal.VectorFunction, 1)
	input[0] = signal.NewInput(func(arg1 float64) float64 {
		return 1.
	}, B)
	stateSpaceModel := ssm.NewIntegratorChain(N, stageGain, input)
	Q := mat.NewDense(N, N, nil)
	Q.Outer(sigma_u, B, B)

	var Rinv mat.Dense
	Rinv.Scale(1./sigma_z, gonumExtensions.Eye(N, N, 0))

	Vf := mat.NewDense(N, N, nil)
	for index := 0; index < N; index++ {
		Vf.Set(index, index, 10000)
	}
	var Res mat.Matrix
	var tmp1 mat.Dense
	var R mat.Dense

	tmp1.Mul(&Rinv, stateSpaceModel.C)
	R.Mul(stateSpaceModel.C.T(), &tmp1)
	Res = care(stateSpaceModel.A.T(), &R, Q, Vf, MatrixFactorization{})
	fmt.Println(mat.Formatted(Res))
}
