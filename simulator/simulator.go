package simulator

import (
	"sync"

	"github.com/hammal/ssm"
	"gonum.org/v1/gonum/mat"
)

// The Simulator interface describes a simulator interface
type Simulator interface {
	// Generate a new state space model with an intalState
	New(ssm ssm.CSSM, initalState mat.Vector, Ts float64)

	// Run the simulator with input function until time t.
	Run(inputs []Signal, t float64)

	// Compute input
	compute(input Signal, t_from float64, t_to float64, returnChannel chan mat.VecDense, sync sync.WaitGroup)
}

// Simulator struct (this needs a better name)
type generallSimulator struct {
	state     *mat.VecDense
	ssm       ssm.CSSM
	t, Ts, Fs float64
	Ad        mat.Dense
}

// The initalisation function for zero order hold simulator
func (sim generallSimulator) New(ssm ssm.CSSM, initalState mat.Vector, Ts float64) {
	sim.ssm = ssm
	sim.state = mat.VecDenseCopyOf(initalState)
	sim.Ts = Ts
	sim.Fs = 1 / Ts
	sim.t = 0
}

func (sim *generallSimulator) Run(inputs []Signal, t float64) {
	// This is the return channel for all vector corresponding to inputs
	returnChannel := make(chan mat.VecDense)
	var wg sync.WaitGroup

	for _, input := range inputs {
		wg.Add(1)
		go sim.compute(input, sim.t, t, returnChannel, wg)
	}

	sim.state.MulVec(&sim.Ad, sim.state)

	// Summaraise all contributions into state vector
	wg.Add(1)
	go func() {
		defer wg.Done()
		for res := range returnChannel {
			sim.state.AddVec(sim.state, &res)
		}
	}()

	// Wait until all input contributions have been computed
	wg.Wait()
}

func (sim *generallSimulator) compute(input Signal, t_from float64, t_to float64, returnChannel chan mat.VecDense, wg sync.WaitGroup) {
	defer wg.Done()
	var tempV mat.VecDense
	switch kind := input.kind; kind {
	case "zero order hold":
		// This is Ad b
		tempV.MulVec(&sim.Ad, input.vector)
		// Add the input evaluated at t_from
		tempV.ScaleVec(input.signal(t_from), &tempV)
		// To be implemented
	default:
		var numberOfIntegrationSteps int = 1000
		// Do the Runge Kutta method of 4th order
		sim.rk4(*mat.NewVecDense(sim.ssm.Order, nil), t_from, t_to, numberOfIntegrationSteps)
	}
	// Return the resulting vector to return channel
	returnChannel <- tempV
}

func (sim generallSimulator) f(BU *mat.VecDense, y *mat.VecDense) {
	y.MulVec(sim.ssm.A, y)
	y.AddVec(y, BU)
}

func (sim generallSimulator) rk4(y0 *mat.VecDense, u *Signal, t0, t1 float64, N int) mat.VecDense {
	res := y0
	// step size
	var h float64 = (t1 - t0) / float64(N)
	t := t0

	for index := 0; index < N; index++ {
		k1 := 1.
	}

	return y0
}
