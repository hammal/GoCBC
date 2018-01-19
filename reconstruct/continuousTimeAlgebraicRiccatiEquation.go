package reconstruct

import (
	"fmt"

	"gonum.org/v1/gonum/mat"
)

type Recursion struct {
	precision  float64
	stepLength float64
}

type EigenDecomposition struct{}

// care Computes the continuous algebraic Riccati equation as shown in [Optimal Solution to Matrix Riccati Equation – For Kalman Filter Implementation](http://cdn.intechopen.com/pdfs/39345/intech-optimal_solution_to_matrix_riccati_equation_for_kalman_filter_implementation.pdf).
//
// X' = A^TX + XA - X B^H Rinv B X + Q
//
// TODO: This care does not work properly!
func care(A, B, Rinv, Q, X mat.Matrix, method interface{}) {
	switch met := method.(type) {
	case Recursion:
		var (
			Jacobian mat.Dense
			tmp1     mat.Dense
			tmp2     mat.Dense
		)
		// A^T X
		Jacobian.Mul(A.T(), X.T())

		// X A
		tmp1.Mul(X, A)
		Jacobian.Add(&Jacobian, &tmp1)

		// X B^T Rinv B X
		tmp2.Mul(B, X)
		tmp2.Mul(Rinv, &tmp2)
		tmp1.Mul(B.T(), &tmp2)
		tmp1.Mul(X, &tmp1)
		Jacobian.Sub(&Jacobian, &tmp1)

		// Q
		Jacobian.Add(&Jacobian, Q)

		// Check norm
		norm := mat.Norm(&Jacobian, 2)

		// Scale with step length
		Jacobian.Scale(met.stepLength, &Jacobian)

		// Add the Jacobian to the current estimate
		xtmp := X.(*mat.Dense)
		xtmp.Add(xtmp, &Jacobian)

		// check if negative determinant
		// if mat.Det(xtmp) < 0 {
		// 	m, _ := xtmp.Dims()
		// 	tmp1.Scale(-1., gonumExtensions.Eye(m, m, 0))
		// 	xtmp.Mul(&tmp1, xtmp)
		// }

		fmt.Println(mat.Formatted(xtmp))

		// Check if derivative has converged
		if norm > met.precision {
			// If not free memory
			Jacobian.Reset()
			tmp1.Reset()
			tmp2.Reset()
			// Reduce step size
			met.stepLength /= 1.
			// Recursively call yourself
			care(A, B, Rinv, Q, X, method)
		}

	case EigenDecomposition:

		n, _ := A.Dims()
		n2, _ := Rinv.Dims()
		tmp0 := mat.NewDense(n, n2, nil)
		tmp00 := mat.NewDense(n, n, nil)
		tmp1 := mat.NewDense(n, 2*n, nil)
		tmp2 := mat.NewDense(n, 2*n, nil)
		psi := mat.NewDense(2*n, 2*n, nil)
		// tmp1 = [A Q]
		tmp1.Augment(Q, A.T())
		tmp1.Scale(-1, tmp1)
		// tmp0 = B^HR^{-1}B
		tmp0.Mul(B.T(), Rinv)
		tmp00.Mul(tmp0, B)
		tmp00.Scale(-1, tmp00)
		// tmp2 = [B^HR^{-1}B -A]
		tmp2.Augment(A, tmp00)
		// psi = [A, Q; B^HR^{-1}B -A]
		psi.Stack(tmp2, tmp1)

		eigen := mat.Eigen{}

		eigen.Factorize(psi, false, true)

		fmt.Println(mat.Formatted(eigen.Vectors()))
		for _, number := range eigen.Values(nil) {
			fmt.Printf("%e, ", real(number))
		}
		fmt.Print("\n")

		// fmt.Printf("%v\n%v\n%v\n", mat.Formatted(tmp0), mat.Formatted(psi), mat.Formatted(tmp1))

		// This represents the system
		// [x_1(t); x_2(t)]' = psi [x_1(t); x_2(t)]
		// which is ordinary differential equation system with the analytical Solution
		// [x_1(t); x_2(t)] = e^(psi t) [x_1(0), x_2(0)]

		// Compute solution for some large time constant t
		// psi = e^(psi t)
		t := 1e1
		psi.Scale(t, psi)
		psi.Exp(psi)
		// fmt.Printf("%v\n%v\n%v\n", mat.Formatted(tmp0), mat.Formatted(psi), mat.Formatted(tmp1))

		// create [I; I] initial vector
		tmp0.Reset()
		tmp1.Reset()

		ones := make([]float64, n)
		for index := range ones {
			ones[index] = 1.
		}
		tmp3 := mat.NewDiagonal(n, ones)
		tmp0.Stack(tmp3, tmp3)
		tmp1.Mul(psi, tmp0)

		// Alternatively this could also be done using the ode library (could have numerically
		// better properties).

		// Compute solution by solving matrix factorization
		tmpX := X.(*mat.Dense)
		tmpX.Solve(tmp1.Slice(0, n, 0, n), tmp1.Slice(n, 2*n, 0, n))

	default:
		meth := Recursion{
			stepLength: 1e-3,
			precision:  1e-12,
		}
		care(A, B, Rinv, Q, X, meth)
	}
}
