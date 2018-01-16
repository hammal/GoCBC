// Package ode is a ordinary differential equation library that implements the
// Runge-Kutta methods https://en.wikipedia.org/wiki/Runge–Kutta_methods.
// The package requires a state space model, from the https://github.com/hammal/stateSpaceModel
// package, to describe the system.
package ode

import (
	"github.com/hammal/adc/ssm"
	"gonum.org/v1/gonum/mat"
)

// RungeKutta holds the butcherTableau which describes the Runge Kutta method.
type RungeKutta struct {
	Description butcherTableau
}

// Compute the update for a Runge-Kutta system based on a current value at t = from
// , a target time t = to, a initial value x(t=from) = value and a system model ode.
// When algorithm is finished the result is copied into the value.
func (rk RungeKutta) Compute(from, to float64, value *mat.VecDense, ode ssm.StateSpaceModel) {
	// State order
	M, _ := value.Dims()
	// The precomputed derivative points
	K := make([]*mat.VecDense, rk.Description.stages)
	// Step length
	h := to - from
	for index := range K {
		// Initalise an intermediate vector
		tempV := mat.NewVecDense(M, nil)
		tempV.CopyVec(value)
		// Compute the relevant vector by combining previously computed derivate points
		// according to Butcher Tableau.
		for index2, a := range rk.Description.rungeKuttaMatrix[index] {
			tempV.AddScaledVec(tempV, h*a, K[index2])
		}
		// Insert the new derivate point
		K[index] = ode.StateDerivative(from+h*rk.Description.nodes[index], tempV)
	}

	// Sum up the different contributions with relevant weights.
	for index, k := range K {
		value.AddScaledVec(value, h*rk.Description.weights[index], k)
	}
}

// NewRK4 function returns a forth order Runge-Kutta object
func NewRK4() *RungeKutta {
	var temp butcherTableau
	temp.stages = 4
	temp.nodes = []float64{0, 0.5, 0.5, 1}
	temp.weights = []float64{1. / 6, 1. / 3, 1. / 3, 1. / 6}
	temp.rungeKuttaMatrix = [][]float64{
		nil,
		{1. / 2},
		{0, 1. / 2},
		{0, 0, 1.},
	}
	rk := RungeKutta{temp}
	return &rk
}

// NewEulerMethod returns a pointer to a Runge-Kutta that does the Euler method.
func NewEulerMethod() *RungeKutta {
	var temp butcherTableau
	temp.stages = 1
	temp.nodes = []float64{0}
	temp.weights = []float64{1}
	rk := RungeKutta{temp}
	return &rk
}

// butcherTableau which describes the approximate solution, see https://en.wikipedia.org/wiki/Runge–Kutta_methods.
type butcherTableau struct {
	stages           int
	weights          []float64
	nodes            []float64
	rungeKuttaMatrix [][]float64
}
