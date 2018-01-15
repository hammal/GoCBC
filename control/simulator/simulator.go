// Package simulator is an interface for simulating different ADC circuits.
// The general idea is given a state space model, see github.com/hammal/ssm, and
// an array of VectorFunctions, these should be seen as real valued function with
// an associated input vector, the simulator solves the associated ordinary differential
// equations and outputs an array of bitstreams...
package simulator

import (
	"sync"

	"github.com/hammal/ode"
	"github.com/hammal/ssm"
	"gonum.org/v1/gonum/mat"
)

// Simulator interface
type Simulator interface {
	GetObservations() [][]float64
	Simulate(t0, t1 float64, state *mat.VecDense, additionalInputs []*ssm.SSM)
}

// simulator type
type simulator struct {
	Ad           mat.Matrix
	observations [][]float64
	index        int
	signalSsms   []*ssm.SSM
}

// GetObservations returns the array of observations
func (sim simulator) GetObservations() [][]float64 {
	return sim.observations
}

func computeFunc(from, to float64, input *ssm.SSM, returnChannel chan mat.VecDense, sync *sync.WaitGroup) {
	defer sync.Done()
	o := ode.NewRK4()
	var tmp mat.VecDense
	o.Compute(from, to, tmp, input)
	returnChannel <- *tmp
}

func (sim *simulator) Simulate(t0, t1 float64, state *mat.VecDense, sm []*ssm.SSM) {
	var wg sync.WaitGroup
	returnChannel := make(chan mat.VecDense)

	wg.Add(len(sim.signalSsms) + len(sm))

	for _, localSSM := range sim.signalSsms {
		go computeFunc(t0, t1, localSSM, returnChannel, &wg)
	}

	for _, controlSSM := range sm {
		go computeFunc(t0, t1, controlSSM, returnChannel, &wg)
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
	N, _ := sim.signalSsms[0].C.Dims()
	tmp := mat.NewVecDense(N, nil)
	tmp.MulVec(sim.signalSsms[0].C, state)
	for i := 0; i < N; i++ {
		sim.observations[sim.index][i] = tmp.AtVec(i)
	}
}

// NewLinearStateSpaceSimulator returns an initialized simulator
func NewLinearStateSpaceSimulator(sys System, sm []ssm.LinearStateSpaceModel) Simulator {

	// Estimate Ad
	r, c := sm[0].A.Dims()
	Ad := mat.NewDense(r, c, nil)
	Ad.Scale(sys.Ts, sm[0].A)
	Ad.Exp(Ad)

	// Initialize observation matrix
	r, _ = sm[0].C.Dims()
	obs := make([][]float64, sys.N)
	for i := int64(0); i < sys.N; i++ {
		obs[i] = make([]float64, r)
	}

	// Create type with index 0
	return &simulator{Ad, obs, 0, sm}
}
