package samplingnetwork

import (
	"fmt"
	"math"

	"gonum.org/v1/gonum/mat"
)

// Here we present some shorthand functions to generate common topologies

// IntegratorBlock is the fundamental integrator block.
func IntegratorBlock(gain float64) SamplingNetwork {
	system := LinearSystem{
		A: mat.NewDense(1, 1, nil),
		B: mat.NewVecDense(1, []float64{gain}),
		C: mat.NewVecDense(1, []float64{1}),
	}
	control := make([]Control, 1)
	control[0] = &AnalogSwitch{
		vector: mat.NewVecDense(1, []float64{-gain}),
	}
	return SamplingNetwork{
		System:  system,
		Control: control,
	}
}

// OscillatorBlock is the fundamental oscillator where we have made the two states
// equally amplified and the control goes into each state. We call this a symmetric
// oscillator block.
func OscillatorBlock(gain, resonanceFrequency float64) SamplingNetwork {
	system := LinearSystem{
		A: mat.NewDense(2, 2, []float64{0, -2. * math.Pi * resonanceFrequency, 2. * math.Pi * resonanceFrequency, 0}),
		// A: mat.NewDense(2, 2, []float64{0, -2. * math.Pi * resonanceFrequency * 2. * math.Pi * resonanceFrequency, 1, 0}),
		// A: mat.NewDense(2, 2, []float64{0, 1, -2. * math.Pi * resonanceFrequency * 2. * math.Pi * resonanceFrequency, 0}),
		B: mat.NewDense(2, 2, []float64{gain, 0, 0, gain}),
		C: mat.NewDense(2, 2, []float64{1, 0, 0, 1}),
	}
	control := make([]Control, 2)

	control[0] = &Oscillator{
		frequency: resonanceFrequency,
		phase:     0.,
		vector:    mat.NewVecDense(2, []float64{gain, 0}),
	}

	control[1] = &Oscillator{
		frequency: resonanceFrequency,
		phase:     math.Pi / 2.,
		vector:    mat.NewVecDense(2, []float64{gain, 0}),
	}

	return SamplingNetwork{
		System:  system,
		Control: control,
	}
}

func SplitBlock(systems []SamplingNetwork) SamplingNetwork {
	if len(systems) == 2 {

		// Check if splitting is possible (same number of inputs)
		if systems[0].System.InputSpaceOrder() != systems[1].System.InputSpaceOrder() {
			panic("The input dimensions must be equal in order to create split block")
		}

		// Update LinearSystem
		systemOrder0 := systems[0].System.StateSpaceOrder()
		systemOrder1 := systems[1].System.StateSpaceOrder()

		// New A matrix
		var tmpU, tmpL, A mat.Dense
		tmpU.Augment(systems[0].System.A, mat.NewDense(systemOrder0, systemOrder1, nil))
		tmpL.Augment(mat.NewDense(systemOrder1, systemOrder0, nil), systems[1].System.A)
		A.Stack(&tmpU, &tmpL)

		tmpU.Reset()
		tmpL.Reset()

		// New C Matrix
		var C mat.Dense
		tmpU.Augment(systems[0].System.C, mat.NewDense(systems[0].System.OutputSpaceOrder(), systemOrder1, nil))
		tmpL.Augment(mat.NewDense(systems[1].System.OutputSpaceOrder(), systemOrder0, nil), systems[1].System.C)
		C.Stack(&tmpU, &tmpL)

		// New B Matrix
		var B mat.Dense
		B.Stack(systems[0].System.B, systems[1].System.B)

		// Update Controls
		var controls []Control

		for index := range systems[0].Control {
			ctrlVec := systems[0].Control[index].GetVector()
			tmpVec := mat.NewVecDense(systemOrder0+systemOrder1, nil)
			for row := 0; row < systemOrder0; row++ {
				tmpVec.SetVec(row, ctrlVec.AtVec(row))
			}
			systems[0].Control[index].SetVector(tmpVec)
			controls = append(controls, systems[0].Control[index])
		}

		for index := range systems[1].Control {
			ctrlVec := systems[1].Control[index].GetVector()
			tmpVec := mat.NewVecDense(systemOrder0+systemOrder1, nil)
			for row := systemOrder0; row < systemOrder0+systemOrder1; row++ {
				tmpVec.SetVec(row, ctrlVec.AtVec(row-systemOrder0))
			}
			systems[1].Control[index].SetVector(tmpVec)
			controls = append(controls, systems[1].Control[index])
		}

		return SamplingNetwork{
			System: LinearSystem{
				A: &A,
				B: &B,
				C: &C,
			},
			Control: controls,
		}

	}
	// If not 2 systems left pop the first and split with the recursive command
	tmpSysArray := make([]SamplingNetwork, 2)
	tmpSysArray[1] = SplitBlock(systems[1:])
	tmpSysArray[0] = systems[0]
	return SplitBlock(tmpSysArray)

}

func MergeBlock(systems []SamplingNetwork) SamplingNetwork {
	if len(systems) == 2 {

		// Check if Merge is possible
		if systems[0].System.OutputSpaceOrder() != systems[1].System.OutputSpaceOrder() {
			panic("The output dimensions must be equal in order to create a merge block")
		}

		// Update LinearSystem
		systemOrder0 := systems[0].System.StateSpaceOrder()
		systemOrder1 := systems[1].System.StateSpaceOrder()

		// New A matrix
		var tmpU, tmpL, A mat.Dense
		tmpU.Augment(systems[0].System.A, mat.NewDense(systemOrder0, systemOrder1, nil))
		tmpL.Augment(mat.NewDense(systemOrder1, systemOrder0, nil), systems[1].System.A)
		A.Stack(&tmpU, &tmpL)

		tmpU.Reset()
		tmpL.Reset()

		// New B Matrix
		var B mat.Dense
		tmpU.Augment(systems[0].System.B, mat.NewDense(systemOrder0, systems[1].System.InputSpaceOrder(), nil))
		tmpL.Augment(mat.NewDense(systemOrder1, systems[0].System.InputSpaceOrder(), nil), systems[1].System.B)
		B.Stack(&tmpU, &tmpL)

		// New C Matrix
		var C mat.Dense
		C.Augment(systems[0].System.C, systems[1].System.C)

		// Update Controls
		var controls []Control

		for index := range systems[0].Control {
			ctrlVec := systems[0].Control[index].GetVector()
			tmpVec := mat.NewVecDense(systemOrder0+systemOrder1, nil)
			for row := 0; row < systemOrder0; row++ {
				tmpVec.SetVec(row, ctrlVec.AtVec(row))
			}
			systems[0].Control[index].SetVector(tmpVec)
			controls = append(controls, systems[0].Control[index])
		}

		for index := range systems[1].Control {
			ctrlVec := systems[1].Control[index].GetVector()
			tmpVec := mat.NewVecDense(systemOrder0+systemOrder1, nil)
			for row := systemOrder0; row < systemOrder0+systemOrder1; row++ {
				tmpVec.SetVec(row, ctrlVec.AtVec(row-systemOrder0))
			}
			systems[1].Control[index].SetVector(tmpVec)
			controls = append(controls, systems[1].Control[index])
		}

		return SamplingNetwork{
			System: LinearSystem{
				A: &A,
				B: &B,
				C: &C,
			},
			Control: controls,
		}

	}
	// If not 2 systems left pop the first and split with the recursive command
	tmpSysArray := make([]SamplingNetwork, 2)
	tmpSysArray[1] = MergeBlock(systems[1:])
	tmpSysArray[0] = systems[0]
	return MergeBlock(tmpSysArray)

}

func ParallelBlock(systems []SamplingNetwork) SamplingNetwork {
	if len(systems) == 2 {

		// No restrictions on output and input

		// Update LinearSystem
		systemOrder0 := systems[0].System.StateSpaceOrder()
		systemOrder1 := systems[1].System.StateSpaceOrder()

		// New A matrix
		var tmpU, tmpL, A mat.Dense
		tmpU.Augment(systems[0].System.A, mat.NewDense(systemOrder0, systemOrder1, nil))
		tmpL.Augment(mat.NewDense(systemOrder1, systemOrder0, nil), systems[1].System.A)
		A.Stack(&tmpU, &tmpL)

		tmpU.Reset()
		tmpL.Reset()

		// New B Matrix
		var B mat.Dense
		tmpU.Augment(systems[0].System.B, mat.NewDense(systemOrder0, systems[1].System.InputSpaceOrder(), nil))
		tmpL.Augment(mat.NewDense(systemOrder1, systems[0].System.InputSpaceOrder(), nil), systems[1].System.B)
		B.Stack(&tmpU, &tmpL)

		// New C Matrix
		var C mat.Dense
		tmpU.Augment(systems[0].System.C, mat.NewDense(systems[0].System.OutputSpaceOrder(), systemOrder1, nil))
		tmpL.Augment(mat.NewDense(systems[1].System.OutputSpaceOrder(), systemOrder0, nil), systems[1].System.C)
		C.Stack(&tmpU, &tmpL)

		// Update Controls
		var controls []Control

		for index := range systems[0].Control {
			ctrlVec := systems[0].Control[index].GetVector()
			tmpVec := mat.NewVecDense(systemOrder0+systemOrder1, nil)
			for row := 0; row < systemOrder0; row++ {
				tmpVec.SetVec(row, ctrlVec.AtVec(row))
			}
			systems[0].Control[index].SetVector(tmpVec)
			controls = append(controls, systems[0].Control[index])
		}

		for index := range systems[1].Control {
			ctrlVec := systems[1].Control[index].GetVector()
			tmpVec := mat.NewVecDense(systemOrder0+systemOrder1, nil)
			for row := systemOrder0; row < systemOrder0+systemOrder1; row++ {
				tmpVec.SetVec(row, ctrlVec.AtVec(row-systemOrder0))
			}
			systems[1].Control[index].SetVector(tmpVec)
			controls = append(controls, systems[1].Control[index])
		}

		return SamplingNetwork{
			System: LinearSystem{
				A: &A,
				B: &B,
				C: &C,
			},
			Control: controls,
		}

	}
	// If not 2 systems left pop the first and split with the recursive command
	tmpSysArray := make([]SamplingNetwork, 2)
	tmpSysArray[1] = ParallelBlock(systems[1:])
	tmpSysArray[0] = systems[0]
	return ParallelBlock(tmpSysArray)

}

func SeriesBlock(systems []SamplingNetwork) SamplingNetwork {
	if len(systems) == 2 {

		// Check if series connection is possible.
		if systems[0].System.OutputSpaceOrder() != systems[1].System.InputSpaceOrder() {
			panic("The first systems output space must be equal to the second systems input space.")
		}

		// Update LinearSystem
		systemOrder0 := systems[0].System.StateSpaceOrder()
		systemOrder1 := systems[1].System.StateSpaceOrder()

		// New A matrix
		var tmpU, tmpL, A, CB mat.Dense
		tmpU.Augment(systems[0].System.A, mat.NewDense(systemOrder0, systemOrder1, nil))
		CB.Mul(systems[1].System.B, systems[0].System.C)
		tmpL.Augment(&CB, systems[1].System.A)
		A.Stack(&tmpU, &tmpL)

		tmpU.Reset()
		tmpL.Reset()

		// New B Matrix
		var B mat.Dense
		B.Stack(systems[0].System.B, mat.NewDense(systemOrder1, systems[0].System.InputSpaceOrder(), nil))

		// New C Matrix
		var C mat.Dense
		C.Augment(mat.NewDense(systems[1].System.OutputSpaceOrder(), systemOrder0, nil), systems[1].System.C)

		// Update Controls
		var controls []Control

		for index := range systems[0].Control {
			ctrlVec := systems[0].Control[index].GetVector()
			tmpVec := mat.NewVecDense(systemOrder0+systemOrder1, nil)
			for row := 0; row < systemOrder0; row++ {
				tmpVec.SetVec(row, ctrlVec.AtVec(row))
			}
			systems[0].Control[index].SetVector(tmpVec)
			controls = append(controls, systems[0].Control[index])
		}

		for index := range systems[1].Control {
			ctrlVec := systems[1].Control[index].GetVector()
			tmpVec := mat.NewVecDense(systemOrder0+systemOrder1, nil)
			for row := systemOrder0; row < systemOrder0+systemOrder1; row++ {
				tmpVec.SetVec(row, ctrlVec.AtVec(row-systemOrder0))
			}
			systems[1].Control[index].SetVector(tmpVec)
			controls = append(controls, systems[1].Control[index])
		}

		return SamplingNetwork{
			System: LinearSystem{
				A: &A,
				B: &B,
				C: &C,
			},
			Control: controls,
		}

	}
	// If not 2 systems left pop the first and split with the recursive command
	tmpSysArray := make([]SamplingNetwork, 2)
	tmpSysArray[1] = SeriesBlock(systems[1:])
	tmpSysArray[0] = systems[0]
	return SeriesBlock(tmpSysArray)

}

// Here we have not thought about how to adjust controls
// Conjecture: Any feedback system of order N can be stabilized using at most
// N Oscillating controls
func FeedbackBlock(system1, system2 SamplingNetwork) SamplingNetwork {
	// Check if feedforward is possible
	if system1.System.OutputSpaceOrder() != system2.System.InputSpaceOrder() {
		panic("The feedforward path between system1 and system2 does not have the same dimensions")
	}

	// Check if splitting is possible (same number of inputs)
	if system1.System.InputSpaceOrder() != system2.System.OutputSpaceOrder() {
		panic("The feedback path between system2 and system2 does not have the same dimensions")
	}

	// Update LinearSystem
	systemOrder0 := system1.System.StateSpaceOrder()
	systemOrder1 := system2.System.StateSpaceOrder()

	// New A matrix
	var tmpU, tmpL, A, CBf, CBb mat.Dense
	CBb.Mul(system1.System.B, system2.System.C)
	// Negative feedback
	CBb.Scale(-1., &CBb)
	tmpU.Augment(system1.System.A, &CBb)
	CBf.Mul(system2.System.B, system1.System.C)
	tmpL.Augment(&CBf, system2.System.A)
	A.Stack(&tmpU, &tmpL)

	tmpU.Reset()
	tmpL.Reset()

	// New B Matrix
	var B mat.Dense
	B.Stack(system1.System.B, mat.NewDense(systemOrder1, system1.System.InputSpaceOrder(), nil))

	// New C Matrix
	var C mat.Dense
	C.Augment(system1.System.C, mat.NewDense(system1.System.OutputSpaceOrder(), systemOrder1, nil))
	// C.Augment(mat.NewDense(system2.System.OutputSpaceOrder(), systemOrder0, nil), system2.System.C)

	// Update Controls
	var controls []Control

	for index := range system1.Control {
		ctrlVec := system1.Control[index].GetVector()
		tmpVec := mat.NewVecDense(systemOrder0+systemOrder1, nil)
		for row := 0; row < systemOrder0; row++ {
			tmpVec.SetVec(row, ctrlVec.AtVec(row))
		}
		system1.Control[index].SetVector(tmpVec)
		controls = append(controls, system1.Control[index])
	}

	for index := range system2.Control {
		ctrlVec := system2.Control[index].GetVector()
		tmpVec := mat.NewVecDense(systemOrder0+systemOrder1, nil)
		for row := systemOrder0; row < systemOrder0+systemOrder1; row++ {
			tmpVec.SetVec(row, ctrlVec.AtVec(row-systemOrder0))
		}
		system2.Control[index].SetVector(tmpVec)
		controls = append(controls, system2.Control[index])
	}

	return SamplingNetwork{
		System: LinearSystem{
			A: &A,
			B: &B,
			C: &C,
		},
		Control: controls,
	}
}

func MultiPlexer(system SamplingNetwork, maps mat.Matrix) SamplingNetwork {
	M, _ := maps.Dims()
	_, N := system.System.B.Dims()

	if M != N {
		panic(fmt.Sprintf("Not a valid mapping for this Sampling Network. B = \n%v\nmaps = \n%v\n", mat.Formatted(system.System.B), mat.Formatted(maps)))
	}
	var tmpB mat.Dense
	tmpB.Mul(system.System.B, maps)
	system.System.B = &tmpB

	// Let's not forget that the controls need to be adjusted
	// TODO way of adjusting control weights.
	return system
}

func DeMultiPlexer(system SamplingNetwork, maps mat.Matrix) SamplingNetwork {
	_, N := maps.Dims()
	M, _ := system.System.C.Dims()

	if M != N {
		panic(fmt.Sprintf("Not a valid mapping for this Sampling Network. C = \n%v\nmaps = \n%v\n", mat.Formatted(system.System.B), mat.Formatted(maps)))
	}
	var tmpC mat.Dense
	tmpC.Mul(maps, system.System.C)
	system.System.C = &tmpC
	return system
}
