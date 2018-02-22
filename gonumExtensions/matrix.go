package gonumExtensions

import (
	"math"

	"gonum.org/v1/gonum/mat"
)

// Ones returns a (m by n) matrix filled with ones
func Ones(m, n int) mat.Matrix {
	return Full(m, n, 1.)
}

// Full returns a (m by n) matrix filled with value
func Full(m, n int, value float64) mat.Matrix {
	data := make([]float64, m*n)
	for index := range data {
		data[index] = value
	}
	tmp := mat.NewDense(m, n, data)
	return tmp
}

// Eye returns a band matrix
func Eye(m, n, k int) mat.Matrix {
	if k == 0 {
		data := make([]float64, int(math.Min(float64(m), float64(n))))
		for entry := range data {
			data[entry] = 1
		}
		return mat.NewDiagonalRect(m, n, data)
	}
	tmp := mat.NewDense(m, n, nil)
	for row := 0; row < m; row++ {
		for column := 0; column < n; column++ {
			if column == (row + k) {
				tmp.Set(row, column+k, 1)
			}
		}
	}
	return tmp
	panic("Not yet implemented")
}

// NANORIF checks if there are any NAN or INF in matrix
func NANORINF(matrix mat.Matrix) bool {
	m, n := matrix.Dims()
	for row := 0; row < m; row++ {
		for col := 0; col < n; col++ {
			if math.IsNaN(matrix.At(row, col)) || math.IsInf(matrix.At(row, col), 0) {
				return true
			}
		}
	}
	return false
}

// Vectorize vectorizes a matrix column wise.
func Vectorize(matrix mat.Matrix) mat.Vector {
	M, N := matrix.Dims()
	res := mat.NewVecDense(M*N, nil)
	for column := 0; column < N; column++ {
		for row := 0; row < M; row++ {
			res.SetVec(row+M*column, matrix.At(row, column))
		}
	}
	return res
}
