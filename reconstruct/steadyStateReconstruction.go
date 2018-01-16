package reconstruct

import (
	"sync"

	"github.com/hammal/adc/control"
	"github.com/hammal/adc/ssm"
	"gonum.org/v1/gonum/mat"
)

type steadyStateReconstruction struct {
	// Forward steady state dynamics
	Af mat.Dense
	// Backward steady state dynamics
	Ab mat.Dense
	// Input Matrix
	W mat.Dense
	// Control interface
	control control.Control
	// state space Model
	stateSpaceModel *ssm.LinearStateSpaceModel
}

func (rec steadyStateReconstruction) Reconstruction() [][]float64 {
	forwardStates := make([]*mat.VecDense, rec.control.Length())
	backwardStates := make([]*mat.VecDense, rec.control.Length())

	m, _ := rec.Af.Dims()

	forwardStates[0] = mat.NewVecDense(m, nil)
	backwardStates[rec.control.Length()-1] = mat.NewVecDense(m, nil)

	var wg sync.WaitGroup

	wg.Add(2)

	// Execute forward and Backward pass concurrently
	go rec.forwardMessagePassing(forwardStates, &wg)
	go rec.backwardMessagePassing(backwardStates, &wg)

	// Wait until done
	wg.Wait()

	// Estimate input
	return rec.inputEstimation(forwardStates, backwardStates)
}

func (rec steadyStateReconstruction) forwardMessagePassing(res []*mat.VecDense, wait *sync.WaitGroup) {
	// tell wait group when done
	defer wait.Done()
	// Forward recursion
	for index := 0; index < len(res)-1; index++ {
		res[index+1].MulVec(&rec.Af, res[index])
		for _, ctrl := range rec.control.GetForwardControlFilterContribution(index) {
			res[index+1].AddVec(res[index+1], ctrl)
		}
	}
}

func (rec steadyStateReconstruction) backwardMessagePassing(res []*mat.VecDense, wait *sync.WaitGroup) {
	// tell wait group when done
	defer wait.Done()
	// Backward recursion
	for index := len(res) - 1; index > 0; index-- {
		res[index-1].MulVec(&rec.Ab, res[index])
		for _, ctrl := range rec.control.GetBackwardControlFilterContribution(index) {
			res[index-1].AddVec(res[index-1], ctrl)
		}
	}
}

func (rec steadyStateReconstruction) inputEstimation(fm []*mat.VecDense, bm []*mat.VecDense) [][]float64 {
	m := rec.stateSpaceModel.InputSpaceOrder()
	n := len(fm)

	// Initialize result vector
	res := make([][]float64, m)
	var tmp mat.VecDense
	for index, _ := range res {
		res[index] = make([]float64, n)
		// (mf + mb) Note  that mb is computed with negation
		tmp.AddVec(fm[index], bm[index])
		tmp.MulVec(&rec.W, &tmp)
		for inp := 0; inp < n; inp++ {
			res[index][inp] = tmp.AtVec(inp)
		}
	}
	return res
}
