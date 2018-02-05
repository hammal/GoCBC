package reconstruct

import (
	"fmt"

	"gonum.org/v1/gonum/mat"
)

type Recursion struct {
	precision  float64
	stepLength float64
}

type NewtonMethod struct {
	precision float64
}

type MatrixFactorization struct{}

// care Computes the continuous algebraic Riccati equation as shown in [Optimal Solution to Matrix Riccati Equation – For Kalman Filter Implementation](http://cdn.intechopen.com/pdfs/39345/intech-optimal_solution_to_matrix_riccati_equation_for_kalman_filter_implementation.pdf).
//
// X' = A^TX + XA - X R X + Q
//
func care(A, R, Q, X mat.Matrix, method interface{}) mat.Matrix {
	switch met := method.(type) {
	case Recursion:
		var (
			Jacobian mat.Dense
			tmp1     mat.Dense
			tmp2     mat.Dense
			xtmp     mat.Dense
		)
		// A^T X
		Jacobian.Mul(A.T(), X)

		// X A
		tmp1.Mul(X, A)
		Jacobian.Add(&Jacobian, &tmp1)

		// X B^T Rinv B X
		tmp2.Mul(R, X)
		tmp1.Mul(X, &tmp2)
		Jacobian.Sub(&Jacobian, &tmp1)

		// Q
		Jacobian.Add(&Jacobian, Q)

		// Check norm
		norm := mat.Norm(&Jacobian, 2)

		// Scale with step length
		Jacobian.Scale(met.stepLength, &Jacobian)

		// Add the Jacobian to the current estimate
		xtmp.Add(X, &Jacobian)

		// check if negative determinant
		// if mat.Det(xtmp) < 0 {
		// 	m, _ := xtmp.Dims()
		// 	tmp1.Scale(-1., gonumExtensions.Eye(m, m, 0))
		// 	xtmp.Mul(&tmp1, xtmp)
		// }

		// fmt.Printf("Norm of jacobian = %v and steplength = %v\n", norm, met.stepLength)

		// Check if derivative has converged
		if norm > met.precision && met.stepLength > 0 {
			// If not free memory
			Jacobian.Reset()
			tmp1.Reset()
			tmp2.Reset()
			// Reduce step size
			// Recursively call yourself
			return care(A, R, Q, &xtmp, Recursion{
				precision:  met.precision,
				stepLength: met.stepLength,
			})
		}
		return &xtmp

	case MatrixFactorization:

		n, _ := A.Dims()
		tmp00 := mat.NewDense(n, n, nil)
		tmp1 := mat.NewDense(n, 2*n, nil)
		tmp2 := mat.NewDense(n, 2*n, nil)
		psi := mat.NewDense(2*n, 2*n, nil)

		// tmp1 = [A, -R]
		tmp00.Scale(1., R)
		tmp1.Augment(A, tmp00)

		// tmp2 = [-Q, -A^T]
		var tmpATM mat.Dense
		tmpATM.Scale(-1., A.T())
		tmp2.Augment(Q, &tmpATM)

		// tmpATM.Scale(1., A.T())
		// tmp2.Augment(Q, &tmpATM)
		// tmp2.Scale(-1., tmp2)

		// psi = [A, -R; -Q, -A^T]
		psi.Stack(tmp1, tmp2)
		// fmt.Printf("tmp1 = \n%v\ntmp2 = \n%v\nPsi = \n%v\n", mat.Formatted(tmp1), mat.Formatted(tmp2), mat.Formatted(psi))

		// Eigenvalue decomposition
		eigen := mat.Eigen{}
		eigen.Factorize(psi, false, true)

		// Extract real eigenvalues
		realEigenValues := make([]float64, 2*n)
		for index, value := range eigen.Values(nil) {
			realEigenValues[index] = real(value)
		}
		// fmt.Printf("EigenValues %v\n", eigen.Values(nil))
		// fmt.Printf("Real eigenvalues %v\n", realEigenValues)

		// Sort eigenvalues into positive and negative partion.
		_, swaps := bubblesort(realEigenValues)

		// Create permutation Matrix
		var permutationMatrix mat.Dense
		permutationMatrix.Permutation(2*n, swaps)

		// Get sorted eigenvectors
		var sortedEigenVectors mat.Dense
		sortedEigenVectors.Mul(eigen.Vectors(), &permutationMatrix)

		fmt.Println(mat.Formatted(&permutationMatrix))
		fmt.Println(mat.Formatted(&sortedEigenVectors))

		// U_2 U_1^(-1)
		// where U_2 is the lower half column vector of the Eigenvectors
		// associated with the positive Eigenvalues.
		var tmpX mat.Dense
		U_1 := sortedEigenVectors.Slice(0, n, 0, n)
		U_2 := sortedEigenVectors.Slice(n, 2*n, 0, n)
		tmpX.Inverse(U_1)
		// fmt.Println(mat.Formatted(&tmpX))
		tmpX.Mul(U_2, &tmpX)

		// Check if negative definite
		if mat.Det(&tmpX) < 0 {
			// If so make positive definite
			tmpX.Scale(-1, &tmpX)
		}

		return &tmpX

	case NewtonMethod:
		var (
			RiccatiEquation mat.Dense
			Jacobian        mat.Dense
			tmp1            mat.Dense
			tmp2            mat.Dense
			x               mat.Dense
		)

		norm := 2 * met.precision
		x = *X.(*mat.Dense)

		for norm > met.precision {
			// A^T X
			RiccatiEquation.Mul(A.T(), &x)

			// X A
			tmp1.Mul(&x, A)
			RiccatiEquation.Add(&RiccatiEquation, &tmp1)

			// X R X
			tmp2.Mul(R, &x)
			tmp1.Mul(&x, &tmp2)
			RiccatiEquation.Sub(&RiccatiEquation, &tmp1)

			// Q
			RiccatiEquation.Add(&RiccatiEquation, Q)

			// (A - RX)
			tmp1.Sub(A, &tmp2)
			// X (A - RX)
			tmp2.Mul(&x, &tmp1)
			// (A - RX)^T X
			tmp1.Mul(tmp1.T(), &x)
			// X(A - RX) + (A - RX)^T X
			// Negative Gradient
			Jacobian.Add(&tmp2, &tmp1)
			// Jacobian.Scale(-1, &Jacobian)

			tmp1.Inverse(&Jacobian)
			tmp1.Mul(&tmp1, &RiccatiEquation)

			// Check norm
			norm = mat.Norm(&tmp1, 2)

			// Add the newton method.
			x.Add(&x, &tmp1)

			// Symmetric
			x.Add(&x, x.T())
			x.Scale(0.5, &x)

			fmt.Printf("Norm of Newton update = %v\n", norm)
			fmt.Println(mat.Formatted(&x))

		}
		return &x

	default:
		meth := Recursion{
			stepLength: 1e-3,
			precision:  1e-12,
		}
		return care(A, R, Q, X, meth)
	}
}

func bubblesort(list []float64) ([]float64, []int) {
	var (
		tmpSwap int
		tmpList float64
	)
	// Initiate swapping list as normal
	length := len(list)
	swap := make([]int, length)
	for index := range list {
		swap[index] = index
	}

	swapped := true
	for swapped {
		swapped = false
		for index := 1; index < length; index++ {
			// fmt.Printf("Value low %v and high %v for index %v \n", list[index-1], list[index], index)
			if list[index-1] > list[index] {
				tmpSwap = swap[index]
				tmpList = list[index]

				swap[index] = swap[index-1]
				swap[index-1] = tmpSwap

				list[index] = list[index-1]
				list[index-1] = tmpList

				swapped = true
			}
		}
		length--
	}
	return list, swap
}

//
