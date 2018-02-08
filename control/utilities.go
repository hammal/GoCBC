package control

import (
	"github.com/hammal/adc/ode"
	"gonum.org/v1/gonum/mat"
)

// bitToIndex handy function to convert bit index array into unique index
func bitToIndex(bits []uint) uint {
	var sum uint = 0
	for index := range bits {
		sum += (1 << uint(index)) * bits[index]
	}
	return sum
}

// inverse of bitToIndex
func indexToBits(index uint, length int) []uint {
	bits := make([]uint, length)
	for cont := range bits {
		if ((index >> uint(cont)) & 1) > 0 {
			bits[cont] = 1
		}
	}
	return bits
}

func indexToVec(index uint, length int) mat.Vector {
	bits := indexToBits(index, length)
	res := make([]float64, length)
	for bit := range bits {
		res[bit] = float64(bits[bit]*2 - 1)
	}
	return mat.NewVecDense(length, res)
}

type lazyCache struct {
	cache    []mat.Vector
	computed []bool
	aSwitch  ControlVector
}

func (lc *lazyCache) GetVector(codeWord uint) mat.Vector {
	// Check if in cache
	if lc.computed[codeWord] {
		// If so return cache
		// fmt.Printf("%v in cache\n", codeWord)
		return lc.cache[codeWord]
	} else {
		// fmt.Printf("%v not in cache\n", codeWord)
		// Compute cache element
		lc.cache[codeWord] = lc.aSwitch.GetVector(codeWord)
		// fmt.Printf("Computed vector \n%v\n", mat.Formatted(lc.cache[codeWord]))
		// Set cache index to true
		lc.computed[codeWord] = true
		// Recursively call yourself
		return lc.GetVector(codeWord)
	}
}

func Solve(system ode.DifferentiableSystem, from, to float64, initalState mat.Matrix) mat.Vector {
	o := ode.NewFehlberg45()
	var res mat.Matrix
	if initalState == nil {
		value := mat.NewDense(system.Order(), 1, nil)
		res, _ = o.Compute(from, to, value, system)
	} else {
		res, _ = o.Compute(from, to, initalState, system)
	}
	res2 := res.(*mat.Dense)
	// fmt.Printf("Derivative at t=3 is \n%v\n", mat.Formatted(system.Derivative(3, mat.NewVecDense(system.Order(), nil))))
	// fmt.Printf("Solution of Solve is \n%v\n", mat.Formatted(res))
	return res2.ColView(0)
}

// // zeroOrderHold solves the problem
// // int_0^T_s e^(A(T_s - t)) dt
// // by solving the initial value problem
// // x'(t) = A x(t) + B
// func zeroOrderHold(A mat.Matrix, t float64) mat.Matrix {
// 	// Determine allowed error
// 	const err float64 = 1e-9
//
// 	fmt.Println(mat.Formatted(A))
//
// 	// Define variables
// 	M, _ := A.Dims()
// 	res := mat.NewDense(M, M, nil)
// 	input := make([]signal.VectorFunction, 1)
// 	od := ode.NewFehlberg45()
//
// 	// For each possible unity input
// 	for column := 0; column < M; column++ {
// 		// Create unit vector
// 		tmp := mat.NewVecDense(M, nil)
// 		tmp.SetVec(column, 1.)
// 		input[0] = signal.NewInput(func(arg1 float64) float64 { return 1. }, tmp)
// 		// Assign temporary state space model and initial state vector
// 		tempSSM := ssm.NewLinearStateSpaceModel(A, A, input)
// 		tempState := mat.NewDense(M, 1, nil)
// 		// Solve initial value problem using adaptive Runge-Kutta Fehlberg4(5) method.
// 		// od.AdaptiveCompute(0., t, err, tempState, tempSSM)
// 		tempRes, _ := od.Compute(0., t, tempState, tempSSM)
//
// 		// fmt.Println(mat.Formatted(tempRes))
// 		// fmt.Println(M)
//
// 		// Fill up result matrix for given unit vector.
// 		for row := 0; row < M; row++ {
// 			res.Set(row, column, tempRes.At(row, 0))
// 		}
// 	}
// 	return res
// }
