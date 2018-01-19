package reconstruct

import (
	"fmt"
	"sync"

	"github.com/hammal/adc/control"
	"github.com/hammal/adc/ssm"
	"gonum.org/v1/gonum/mat"
)

// steadyStateReconstruction holds a reconstruction struct that
// solves the message passing based on the idea of that the associated
// covariance matrices converge to a steady state.
type steadyStateReconstruction struct {
	// Forward steady state dynamics
	Af mat.Dense
	// Backward steady state dynamics
	Ab mat.Dense
	// Input Matrix
	W mat.Dense
	// Control interface
	control control.Control
}

// Reconstruction returns the reconstructed estimates based on the associated
// control interface.
//
// The returned data structure is [number of time indices][number of estimates]*float64
func (rec steadyStateReconstruction) Reconstruction() [][]*float64 {

	// Sync group for go routines
	var wg sync.WaitGroup

	// Number of estimates to produce
	n := rec.control.GetLength()

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

	// Initialize result vector
	res := make([][]*float64, n)
	numberOfEstimates, _ := rec.W.Dims()
	for index := range res {
		res[index] = make([]*float64, numberOfEstimates)
	}

	// Keep track of what stages have been computed and commence input estimations if possible
	//
	// These are the two slices for keeping track.
	forwardComputedStates := make([]bool, n)
	backwardComputedStates := make([]bool, n)

	var remainingNumberOfEstimates = n

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

				fmt.Printf("Started ank input estimation from %s for index %v\n", "forward", f)
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

// forwardMessagePassing is a utility function that performs the forward message passing
// and updates the res vector
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

// backwardMessagePassing is a utility function that performs the backward message passing
// and updates the corresponding res vector.
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

// inputEstimate is a utility function that performs the last merging input estimation
// by combining a forward and backward message into estimate.
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

// NewSteadyStateReconstructor returns a Steady-state reconstructor based on the
// control.
func NewSteadyStateReconstructor(cont control.Control, measurementNoiseCovariance, inputNoiseCovariance mat.Matrix, linearStateSpaceModel ssm.LinearStateSpaceModel) *steadyStateReconstruction {
	// This function needs to do the following
	// - Compute steady state covariance matrices
	// 	- Compute filter state dynamics, forward and backward.
	// 		- Call controls PreComputeFilterContributions with state dynamics.
	// 		- Compute the W
	// 		- Compute the Af
	// 		- Compute the Ab
	//

	// variables
	var (
		inverseMeasurementNoiseCovariance *mat.Dense
		order                             int
		ForwardStateDynamics              *mat.Dense
		BackwardStateDynamics             *mat.Dense
		tmpMatrix1                        *mat.Dense
		tmpMatrix2                        *mat.Dense
		W, Af, Ab                         *mat.Dense
		Vf, Vb                            *mat.Dense
		rec                               steadyStateReconstruction
	)

	// Compute state order
	order, _ = linearStateSpaceModel.A.Dims()

	// Compute inverse measurement noise covariance
	inverseMeasurementNoiseCovariance.Inverse(measurementNoiseCovariance)

	// Solve forward and backward steady state covariance function
	Vf = mat.NewDense(order, order, nil)
	careOption := Recursion{
		precision:  1e-9,
		stepLength: 1e-4,
	}
	care(linearStateSpaceModel.A, linearStateSpaceModel.C, inverseMeasurementNoiseCovariance, inputNoiseCovariance, Vf, careOption)

	// Backward recursion with sign changes
	Vb = mat.NewDense(order, order, nil)
	tmpMatrix2.Scale(-1, inverseMeasurementNoiseCovariance)
	tmpMatrix1.Scale(-1, inputNoiseCovariance)
	care(linearStateSpaceModel.A, linearStateSpaceModel.C, tmpMatrix1, tmpMatrix2, Vb, careOption)

	// Reset sizes and clear memory
	tmpMatrix1.Reset()
	tmpMatrix2.Reset()

	// Compute state dynamics
	// Forward: (A - Vf C Sigma_z^(-1) C^T )
	// Backward: -(A + Vb C Sigma_z^(-1) C^T )
	tmpMatrix2.Mul(inverseMeasurementNoiseCovariance, linearStateSpaceModel.C.T())
	tmpMatrix1.Mul(linearStateSpaceModel.C, tmpMatrix2)
	tmpMatrix2.Reset()
	tmpMatrix2.Mul(Vf, tmpMatrix1)
	ForwardStateDynamics.Sub(linearStateSpaceModel.A, tmpMatrix2)

	tmpMatrix2.Reset()
	tmpMatrix2.Mul(Vb, tmpMatrix1)
	tmpMatrix2.Add(linearStateSpaceModel.A, tmpMatrix2)
	BackwardStateDynamics.Scale(-1, tmpMatrix2)

	// Let control initialize the filter contributions
	cont.PreComputeFilterContributions(ForwardStateDynamics, BackwardStateDynamics)

	// Compute input estimate weights
	tmpMatrix1.Reset()
	tmpMatrix1.Add(Vf, Vb)

	tmpMatrix2.Reset()
	tmpMatrix2.Grow(order, len(linearStateSpaceModel.Input))
	for index, input := range linearStateSpaceModel.Input {
		tmpMatrix2.SetCol(index, input.B.RawVector().Data)
	}
	W.Solve(tmpMatrix1, tmpMatrix2)

	// Compute Af
	Af.Scale(cont.GetTs(), ForwardStateDynamics)
	Af.Exp(Af)

	// Compute Bf
	Ab.Scale(cont.GetTs(), BackwardStateDynamics)
	Ab.Exp(Ab)

	// Initialize steady state reconstruction instance
	rec = steadyStateReconstruction{
		*Af,
		*Ab,
		*W,
		cont,
	}

	return &rec
}
