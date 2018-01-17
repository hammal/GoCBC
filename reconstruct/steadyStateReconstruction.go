package reconstruct

import (
	"fmt"
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

func (rec steadyStateReconstruction) Reconstruction() [][]*float64 {

	// Sync group for go routines
	var wg sync.WaitGroup

	// Number of estimates to produce
	n := rec.control.Length()

	// Initialize the state
	forwardStates := make([]*mat.VecDense, n)
	backwardStates := make([]*mat.VecDense, n)

	// Number of states
	m, _ := rec.Af.Dims()

	// Create two channels for reporting back computed states
	forwardReport := make(chan int, 10)
	backwardReport := make(chan int, 10)

	// Initialize the first respectively last state
	forwardStates[0] = mat.NewVecDense(m, nil)
	backwardStates[n-1] = mat.NewVecDense(m, nil)
	// Report that these states have been computed
	forwardReport <- 0
	backwardReport <- n - 1

	// Add forward and backward recursion to wait group
	wg.Add(2)

	// Execute forward and Backward pass concurrently
	go rec.forwardMessagePassing(forwardStates, &wg, forwardReport)
	go rec.backwardMessagePassing(backwardStates, &wg, backwardReport)

	// Number of inputs to be estimated
	numberOfInputs := rec.stateSpaceModel.InputSpaceOrder()

	// Initialize result vector
	res := make([][]*float64, n)
	for index := range res {
		res[index] = make([]*float64, numberOfInputs)
	}

	// Keep track of what stages have been computed and commence input estimations if possible
	//
	// These are the two slices for keeping track.
	forwardComputedStates := make([]bool, n)
	backwardComputedStates := make([]bool, n)

	var remainingNumberOfEstimates int = n

	// Add the number of go routines required to complete input estimation
	wg.Add(remainingNumberOfEstimates)

	// Loop until there are no more estimates to be initiated
	for remainingNumberOfEstimates > 0 {
		// Check if there are any received reports
		select {
		// If there is a forward computation that finished
		case f := <-forwardReport:
			// if state has already been computed for the backward recursion?
			if backwardComputedStates[f] {
				// Start input estimation
				go rec.inputEstimate(res[f], forwardStates[f], backwardStates[f], &wg)
				// Reduce the number of remaining input estimate computations
				remainingNumberOfEstimates--

				fmt.Printf("Started an input estimation from %s for index %v\n", "forward", f)
			} else {
				// Register that this state has been computed
				forwardComputedStates[f] = true
			}
			// If there is a backward computation that finished
		case b := <-backwardReport:
			// check if state has already been computed by the forward recursion?
			if forwardComputedStates[b] {
				// Start input estimation
				go rec.inputEstimate(res[b], forwardStates[b], backwardStates[b], &wg)
				// Reduce the number of remaining estimates
				remainingNumberOfEstimates--

				fmt.Printf("Started an input estimation from %s for index %v\n", "backward", b)
			} else {
				// Register the state as computed
				backwardComputedStates[b] = true
			}
		default:
			fmt.Println("Nothing ready for computing")
		}
	}

	// Wait until all computations done
	wg.Wait()
	return res
}

func (rec steadyStateReconstruction) forwardMessagePassing(res []*mat.VecDense, wait *sync.WaitGroup, report chan<- int) {
	// upon completion tell wait group that you are done and close the report channel
	defer close(report)
	defer wait.Done()
	// Forward recursion
	for index := 0; index < len(res)-1; index++ {
		// Increment the state
		res[index+1].MulVec(&rec.Af, res[index])
		for _, ctrl := range rec.control.GetForwardControlFilterContribution(index) {
			// Check that vector interface is a mat.VecDense type
			c, ok := ctrl.(*mat.VecDense)
			// If not raise panic
			if ok {
				str := fmt.Sprintf("lacking implementation for type %T", ctrl)
				panic(str)
			}
			// Add control contribution to current state
			res[index+1].AddVec(res[index+1], c)
			// Report back that state has been computed
			report <- index + 1
		}
	}
}

func (rec steadyStateReconstruction) backwardMessagePassing(res []*mat.VecDense, wait *sync.WaitGroup, report chan<- int) {
	// upon completion tell wait group that you are done and close the report channel
	defer wait.Done()
	defer close(report)

	// Backward recursion
	for index := len(res) - 1; index > 0; index-- {
		// Increment the state
		res[index-1].MulVec(&rec.Ab, res[index])
		for _, ctrl := range rec.control.GetBackwardControlFilterContribution(index) {
			// Check that vector interface is a mat.VecDense type
			c, ok := ctrl.(*mat.VecDense)
			// If not raise panic
			if ok {
				str := fmt.Sprintf("lacking implementation for type %T", ctrl)
				panic(str)
			}
			// Add control contribution to current state
			res[index-1].AddVec(res[index-1], c)
			// Report back that state has been computed
			report <- index - 1
		}
	}
}

func (rec steadyStateReconstruction) inputEstimate(res []*float64, fm, bm *mat.VecDense, sync *sync.WaitGroup) {
	// Report when function is done
	defer sync.Done()
	// constants
	var (
		tmp      mat.VecDense
		nrInputs int
		tmpRes   float64
	)
	nrInputs = len(res)

	// tmp = (fm - bm)
	tmp.SubVec(fm, bm)
	// tmp = W tmp
	tmp.MulVec(&rec.W, &tmp)
	// Fill up result vector
	for inp := 0; inp < nrInputs; inp++ {
		tmpRes = tmp.AtVec(inp)
		res[inp] = &tmpRes
	}
}
