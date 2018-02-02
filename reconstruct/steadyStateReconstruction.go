package reconstruct

import (
	// "fmt"
	"fmt"
	"sync"

	"github.com/hammal/adc/control"
	"github.com/hammal/adc/gonumExtensions"
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
	// estimate
	estimate [][]float64
	// ForwardMessage
	forwardMessage []mat.VecDense
	// BackwardMessage
	backwardMessage []mat.VecDense
}

// Reconstruction returns the reconstructed estimates based on the associated
// control interface.
//
// The returned data structure is [number of time indices][number of estimates]*float64
func (rec *steadyStateReconstruction) Reconstruction() [][]float64 {

	// Sync group for go routines
	var wg sync.WaitGroup

	// Number of estimates to produce
	n := rec.control.GetLength()

	// Initialize the state
	rec.forwardMessage = make([]mat.VecDense, n)
	rec.backwardMessage = make([]mat.VecDense, n)

	// Number of states
	m, _ := rec.Af.Dims()

	// Create two channels for reporting back computed states
	forwardReport := make(chan int, 100)
	backwardReport := make(chan int, 100)

	// Initialize the first respectively last state
	rec.forwardMessage[0] = *mat.NewVecDense(m, nil)
	rec.backwardMessage[n-1] = *mat.NewVecDense(m, nil)
	// Report that these states have been computed
	forwardReport <- 0
	backwardReport <- n - 1

	// Add forward and backward recursion to wait group
	wg.Add(2)

	// Execute forward and Backward pass concurrently
	go rec.forwardMessagePassing(&wg, forwardReport)
	go rec.backwardMessagePassing(&wg, backwardReport)

	// Keep track of what stages have been computed and commence input estimations if possible
	//
	// These are the two slices for keeping track.
	forwardComputedStates := SafeLedger{ledger: make([]bool, n)}
	backwardComputedStates := SafeLedger{ledger: make([]bool, n)}
	remainingNumberOfEstimates := SafeCounter{v: n}

	// Add the number of go routines required to complete input estimation
	wg.Add(remainingNumberOfEstimates.Get())

	rec.estimate = make([][]float64, n)
	_, nrInputs := rec.W.Dims()
	for tmpEst := range rec.estimate {
		rec.estimate[tmpEst] = make([]float64, nrInputs)
	}

	// Loop until there are no more estimates to be initiated
	for remainingNumberOfEstimates.Get() > 0 {
		// Check if there are any received reports
		select {
		// If there is a forward computation that finished
		case f := <-forwardReport:
			// if state has already been computed for the backward recursion?
			if backwardComputedStates.Get(f) {
				// Start input estimation
				go rec.inputEstimate(f, &wg)
				// Reduce the number of remaining input estimate computations
				remainingNumberOfEstimates.Dec(1)

				// fmt.Printf("Started an input estimation from %s for index %v\n", "forward", f)
			} else {
				// Register that this state has been computed
				forwardComputedStates.Set(f)
			}
			// If there is a backward computation that finished
		case b := <-backwardReport:
			// check if state has already been computed by the forward recursion?
			// fmt.Print(forwardComputedStates[b])
			if forwardComputedStates.Get(b) {
				// Start input estimation
				go rec.inputEstimate(b, &wg)
				// Reduce the number of remaining estimates
				remainingNumberOfEstimates.Dec(1)

				// fmt.Printf("Started an input estimation from %s for index %v\n", "backward", b)
			} else {
				// Register the state as computed
				backwardComputedStates.Set(b)
			}

			// fmt.Println("After switch statement")
		}
	}

	// Wait until all computations done
	wg.Wait()
	return rec.estimate
}

// forwardMessagePassing is a utility function that performs the forward message passing
// and updates the res vector
func (rec *steadyStateReconstruction) forwardMessagePassing(wait *sync.WaitGroup, report chan<- int) {
	// upon completion tell wait group that you are done and close the report channel
	// defer close(report)
	defer wait.Done()
	// Forward recursion
	for index := 0; index < rec.control.GetLength()-1; index++ {
		// Increment the state
		rec.forwardMessage[index+1].MulVec(&rec.Af, &rec.forwardMessage[index])
		ctrl, err := rec.control.GetBackwardControlFilterContribution(index)
		if err != nil {
			panic(err)
		}
		rec.forwardMessage[index+1].AddVec(&rec.forwardMessage[index+1], ctrl)
		report <- index + 1
	}
}

// backwardMessagePassing is a utility function that performs the backward message passing
// and updates the corresponding res vector.
func (rec *steadyStateReconstruction) backwardMessagePassing(wait *sync.WaitGroup, report chan<- int) {
	// upon completion tell wait group that you are done and close the report channel
	defer wait.Done()
	// defer close(report)

	// Backward recursion
	for index := rec.control.GetLength() - 1; index > 0; index-- {
		// Increment the state
		rec.backwardMessage[index-1].MulVec(&rec.Ab, &rec.backwardMessage[index])
		ctrl, err := rec.control.GetBackwardControlFilterContribution(index)
		if err != nil {
			panic(err)
		}
		rec.backwardMessage[index-1].AddVec(&rec.backwardMessage[index-1], ctrl)
		report <- index - 1
		// fmt.Printf("BackwardIndex %v\n", index-1)
	}
}

// inputEstimate is a utility function that performs the last merging input estimation
// by combining a forward and backward message into estimate.
func (rec *steadyStateReconstruction) inputEstimate(resindex int, sync *sync.WaitGroup) {
	// Report when function is done
	defer sync.Done()
	// constants
	var (
		tmp, res mat.VecDense
		nrInputs int
	)
	_, nrInputs = rec.W.Dims()

	// tmp = (fm - bm)
	tmp.SubVec(&rec.forwardMessage[resindex], &rec.backwardMessage[resindex])
	// tmp.CopyVec(&rec.backwardMessage[resindex])
	// tmp.CopyVec(&rec.forwardMessage[resindex])
	// tmp = W tmp
	// fmt.Printf("W^T = \n%v\ntmp = \n%v\n", mat.Formatted(rec.W.T()), mat.Formatted(&tmp))
	res.MulVec(rec.W.T(), &tmp)
	// Fill up result vector
	for inp := 0; inp < nrInputs; inp++ {
		rec.estimate[resindex][inp] = res.AtVec(inp)
	}
	// fmt.Printf("Computed Input estimate at index %v estimate = %v\n", resindex, rec.estimate[resindex])
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
		inverseMeasurementNoiseCovariance mat.Dense
		order                             int
		ForwardStateDynamics              mat.Dense
		BackwardStateDynamics             mat.Dense
		tmpMatrix1                        mat.Dense
		tmpMatrix2                        mat.Dense
		W, Af, Ab                         mat.Dense
		Vf, Vb                            mat.Matrix
		R                                 mat.Dense
		rec                               steadyStateReconstruction
	)

	// Compute state order
	order = linearStateSpaceModel.StateSpaceOrder()

	// Compute inverse measurement noise covariance
	inverseMeasurementNoiseCovariance.Inverse(measurementNoiseCovariance)
	// Solve forward and backward steady state covariance function
	Vf = gonumExtensions.Eye(order, order, 0)
	// careOption := Recursion{
	// 	precision:  1e-1,
	// 	stepLength: 5e-6,
	// }

	R.Mul(&inverseMeasurementNoiseCovariance, linearStateSpaceModel.C.T())
	R.Mul(linearStateSpaceModel.C, &R)

	Vf = care(linearStateSpaceModel.A.T(), &R, inputNoiseCovariance, Vf, MatrixFactorization{})
	// Backward recursion with sign changes
	Vb = gonumExtensions.Eye(order, order, 0)
	tmpMatrix1.Scale(-1, linearStateSpaceModel.A.T())

	Vb = care(&tmpMatrix1, &R, inputNoiseCovariance, Vb, MatrixFactorization{})

	fmt.Printf("Solution To ARE is:\nVf =\n%v\nVb\n%v\n", mat.Formatted(Vf), mat.Formatted(Vb))

	// Compute state dynamics
	// Forward: (A - Vf C Sigma_z^(-1) C^T )
	// Backward: -(A + Vb C Sigma_z^(-1) C^T )
	tmpMatrix1.Mul(Vf, &R)
	tmpMatrix2.Mul(Vb, &R)

	ForwardStateDynamics.Sub(linearStateSpaceModel.A, &tmpMatrix1)

	tmpMatrix2.Add(linearStateSpaceModel.A, &tmpMatrix2)
	BackwardStateDynamics.Scale(-1, &tmpMatrix2)

	fmt.Printf("Forward and Backward Filter Dynamics are:\nAdf = \n%v\nAdb = \n%v\n", mat.Formatted(&ForwardStateDynamics), mat.Formatted(&BackwardStateDynamics))

	// Let control initialize the filter contributions
	cont.PreComputeFilterContributions(&ForwardStateDynamics, &BackwardStateDynamics)

	// Compute input estimate weights
	tmpMatrix1.Reset()
	tmpMatrix1.Add(Vf, Vb)

	tmpMatrix2.Reset()
	tmpMatrix2 = *mat.NewDense(order, linearStateSpaceModel.InputSpaceOrder(), nil)
	// fmt.Printf("Input order %v and B vector \n%v\n", linearStateSpaceModel.InputSpaceOrder(), mat.Formatted(linearStateSpaceModel.Input[0].B))
	for index, input := range linearStateSpaceModel.Input {
		for row := 0; row < order; row++ {
			tmpMatrix2.Set(row, index, -input.B.AtVec(row))
		}
	}
	W.Solve(&tmpMatrix1, &tmpMatrix2)

	// Compute Af
	Af.Scale(cont.GetTs(), &ForwardStateDynamics)
	Af.Exp(&Af)

	// Compute Bf
	Ab.Scale(cont.GetTs(), &BackwardStateDynamics)
	Ab.Exp(&Ab)

	fmt.Printf("PreComputed Matrix\nAf = \n%v\nAb = \n%v\n", mat.Formatted(&Af), mat.Formatted(&Ab))

	// Initialize steady state reconstruction instance
	rec = steadyStateReconstruction{
		Af:      Af,
		Ab:      Ab,
		W:       W,
		control: cont,
	}

	return &rec
}

type SafeLedger struct {
	ledger []bool
	mux    sync.Mutex
}

func (c *SafeLedger) Get(index int) bool {
	c.mux.Lock()
	defer c.mux.Unlock()
	return c.ledger[index]
}

func (c *SafeLedger) Set(index int) {
	c.mux.Lock()
	defer c.mux.Unlock()
	c.ledger[index] = true
}

// SafeCounter is safe to use concurrently.
type SafeCounter struct {
	v   int
	mux sync.Mutex
}

// Inc increments the counter for the given key.
func (c *SafeCounter) Get() int {
	c.mux.Lock()
	defer c.mux.Unlock()
	return c.v
}

// Value returns the current value of the counter for the given key.
func (c *SafeCounter) Dec(number int) {
	c.mux.Lock()
	// Lock so only one goroutine at a time can access the map c.v.
	defer c.mux.Unlock()
	c.v -= number
	// fmt.Println(c.v)
}
