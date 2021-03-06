package ode

import (
	"github.com/hammal/adc/signal"
	"gonum.org/v1/gonum/mat"
)

type NumericalIntegration struct {
	derivative signal.Signal
}

func (ni NumericalIntegration) Order() int { return 1 }

func (nI NumericalIntegration) Derivative(time float64, state mat.Vector) mat.Vector {
	return nI.derivative.Value(time)
}

func (nI NumericalIntegration) Integrate(from, to float64) mat.Matrix {
	_, n := nI.derivative.Value(0).Dims()
	tmpRes := mat.NewVecDense(n, nil)
	var faultTolerance = 1e-9
	o := NewFehlberg45()
	res, err := o.AdaptiveCompute(from, to, faultTolerance, tmpRes, nI)
	if err != nil {
		panic("Adaptive Computation failed")
	}
	return res
}
