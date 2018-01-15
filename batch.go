package adc

import (
	"errors"

	"github.com/hammal/adc/control"
	"github.com/hammal/stateSpaceModel"
)

func NewBatchADC(Ts float64, N int64, stateSpaceModel stateSpaceModel.stateSpaceModel, ct control.Control) *ADC {
	panic(errors.New("Not yet implemented"))
}
