package samplingnetwork

import (
	"math"

	"github.com/hammal/adc/signal"
	"gonum.org/v1/gonum/mat"
)

// Control are responsible of keeping bounded node values in our sampling network.
type Control interface {
	GetState() uint
	GetResponse() signal.VectorFunction
	SetState(float64)
	SetVector(mat.Vector)
	GetVector() mat.Vector
}

// AnalogSwitch is the simplest control mechanism which sends a constant
// signal back to the sampling network.
type AnalogSwitch struct {
	state  uint
	vector mat.Vector
}

func (ctrl AnalogSwitch) GetState() uint {
	return ctrl.state
}

func (ctrl AnalogSwitch) GetResponse() signal.VectorFunction {
	value := (2.*float64(ctrl.state) - 1.)
	var tmpCtrlVec mat.VecDense
	tmpCtrlVec.ScaleVec(value, ctrl.vector)
	return signal.VectorFunction{
		B: &tmpCtrlVec,
		U: func(arg float64) float64 { return 1. },
	}
}

func (ctrl *AnalogSwitch) SetState(value float64) {
	if value > 0 {
		ctrl.state = 1
	} else {
		ctrl.state = 0
	}
}

func (ctrl *AnalogSwitch) SetVector(vec mat.Vector) {
	ctrl.vector = vec
}

func (ctrl AnalogSwitch) GetVector() mat.Vector {
	return ctrl.vector
}

// Oscillators sends an oscillating signal where the frequency and phase offset
// can be set
type Oscillator struct {
	state            uint
	vector           mat.Vector
	frequency, phase float64
}

func (ctrl Oscillator) GetState() uint {
	return ctrl.state
}

func (ctrl Oscillator) GetResponse() signal.VectorFunction {
	value := (2.*float64(ctrl.state) - 1.)
	var tmpCtrlVec mat.VecDense
	tmpCtrlVec.ScaleVec(value, ctrl.vector)
	return signal.VectorFunction{
		B: &tmpCtrlVec,
		U: func(t float64) float64 { return math.Sin(2.*math.Pi*ctrl.frequency*t + ctrl.phase) },
	}
}

func (ctrl *Oscillator) SetState(value float64) {
	if value > 0 {
		ctrl.state = 1
	} else {
		ctrl.state = 0
	}
}

func (ctrl *Oscillator) SetVector(vec mat.Vector) {
	ctrl.vector = vec
}

func (ctrl Oscillator) GetVector() mat.Vector {
	return ctrl.vector
}
