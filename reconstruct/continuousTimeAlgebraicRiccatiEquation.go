package reconstruct

import (
	"gonum.org/v1/gonum/mat"
)

// care Computes the continuous algebraic Riccati equation as shown in [Optimal Solution to Matrix Riccati Equation – For Kalman Filter Implementation](http://cdn.intechopen.com/pdfs/39345/intech-optimal_solution_to_matrix_riccati_equation_for_kalman_filter_implementation.pdf).
//
// X' = AX + XA^H - X B^H Rinv B X + Q
//
// TODO: This care does not work properly!
func care(A, B, Rinv, Q mat.Matrix) mat.Matrix {
	n, _ := A.Dims()
	n2, _ := Rinv.Dims()
	Am := mat.NewDense(n, n, nil)
	tmp0 := mat.NewDense(n, n2, nil)
	tmp00 := mat.NewDense(n, n, nil)
	tmp1 := mat.NewDense(n, 2*n, nil)
	tmp2 := mat.NewDense(n, 2*n, nil)
	psi := mat.NewDense(2*n, 2*n, nil)
	// tmp1 = [A Q]
	tmp1.Augment(A, Q)
	// tmp0 = B^HR^{-1}B
	tmp0.Mul(B.T(), Rinv)
	tmp00.Mul(tmp0, B)
	Am.Scale(-1, A)
	// tmp2 = [B^HR^{-1}B -A]
	tmp2.Augment(tmp00, Am)
	// psi = [A, Q; B^HR^{-1}B -A]
	psi.Stack(tmp1, tmp2)

	// fmt.Printf("%v\n%v\n%v\n", mat.Formatted(tmp0), mat.Formatted(psi), mat.Formatted(tmp1))

	// This represents the system
	// [x_1(t); x_2(t)]' = psi [x_1(t); x_2(t)]
	// which is ordinary differential equation system with the analytical Solution
	// [x_1(t); x_2(t)] = e^(psi t) [x_1(0), x_2(0)]

	// Compute solution for some large time constant t
	// psi = e^(psi t)
	t := 1e2
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
	var X mat.Dense
	X.Solve(tmp1.Slice(0, n, 0, n), tmp1.Slice(n, 2*n, 0, n))

	return &X
}
