// Package simulator is an interface for simulating different ADC circuits.
// The general idea is given a state space model, see github.com/hammal/stateSpaceModel, and
// an array of VectorFunctions, these should be seen as real valued function with
// an associated input vector, the simulator solves the associated ordinary differential
// equations and outputs an array of bitstreams...
package simulate

import (
	"sync"

	"github.com/hammal/adc/ode"
	"github.com/hammal/adc/ssm"
	"gonum.org/v1/gonum/mat"
)

// Simulator interface
type Simulator interface {
	GetObservations() [][]float64
	Simulate(t0, t1 float64, state *mat.VecDense, additionalInputs []*ssm.StateSpaceModel)
}

// simulator type
type simulator struct {
	Ad                     mat.Matrix
	observations           [][]float64
	index                  int
	signalstateSpaceModels []*ssm.StateSpaceModel
}

// GetObservations returns the array of observations
func (sim simulator) GetObservations() [][]float64 {
	return sim.observations
}

func computeFunc(from, to float64, input *ssm.StateSpaceModel, returnChannel chan mat.VecDense, sync *sync.WaitGroup) {
	defer sync.Done()
	o := ode.NewRK4()
	var tmp mat.VecDense
	o.Compute(from, to, &tmp, input)
	returnChannel <- *tmp
}

func (sim *simulator) Simulate(t0, t1 float64, state *mat.VecDense, sm []*ssm.StateSpaceModel) {
	var wg sync.WaitGroup
	returnChannel := make(chan mat.VecDense)

	wg.Add(len(sim.signalstateSpaceModels) + len(sm))

	for _, localstateSpaceModel := range sim.signalstateSpaceModels {
		go computeFunc(t0, t1, localstateSpaceModel, returnChannel, &wg)
	}

	for _, controlstateSpaceModel := range sm {
		go computeFunc(t0, t1, controlstateSpaceModel, returnChannel, &wg)
	}

	go func() {
		wg.Wait()
		close(returnChannel)
	}()

	// Summarizes all contributions into state vector
	state.MulVec(sim.Ad, state)
	for res := range returnChannel {
		state.AddVec(state, &res)
	}

	// Wait until all input contributions have been computed

	// Extract and copy the observations
	N, _ := sim.signalstateSpaceModels[0].C.Dims()
	tmp := mat.NewVecDense(N, nil)
	tmp.MulVec(sim.signalstateSpaceModels[0].C, state)
	for i := 0; i < N; i++ {
		sim.observations[sim.index][i] = tmp.AtVec(i)
	}
}

// NewLinearStateSpaceSimulator returns an initialized simulator
func NewLinearStateSpaceSimulator(order int, Ts float64, sm []ssm.LinearStateSpaceModel) Simulator {

	// Estimate Ad
	r, c := sm[0].A.Dims()
	Ad := mat.NewDense(r, c, nil)
	Ad.Scale(Ts, sm[0].A)
	Ad.Exp(Ad)

	// Initialize observation matrix
	r, _ = sm[0].C.Dims()
	obs := make([][]float64, order)
	for i := int64(0); i < order; i++ {
		obs[i] = make([]float64, r)
	}

	// Create type with index 0
	return &simulator{Ad, obs, 0, sm}
}
