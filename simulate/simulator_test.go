package simulate

import (
	"fmt"
	"testing"

	"github.com/hammal/stateSpaceModel"
	"gonum.org/v1/gonum/mat"
)

// This test does three things
// 1) build input vector
// 2) build state space Model
// 3) test impulse response
func TestNewLinearStateSpaceSimulator(t *testing.T) {
	var n int = 10
	data := make([]float64, n)
	data[0] = 1.
	inp := stateSpaceModel.NewInput(func(arg1 float64) float64 { return 1 }, mat.NewVecDense(n, data))

	sm := make([]stateSpaceModel.LinearStateSpaceModel, 1)
	sm[0] = *stateSpaceModel.NewIntegratorChain(n, 10, inp)
	sys := System{1., 0., 19., 1000, 1, 1}
	newLinearStateSpaceSimulator := NewLinearStateSpaceSimulator(sys, sm)

	state := mat.NewVecDense(n, nil)
	t0, t1 := 0., sys.Ts
	newLinearStateSpaceSimulator.Simulate(t0, t1, state, sm)
	for index := int64(0); index < sys.N; index++ {
		t0 = t1
		t1 += sys.Ts
		newLinearStateSpaceSimulator.Simulate(t0, t1, state, nil)
		fmt.Println(mat.Formatted(state))
	}
}
