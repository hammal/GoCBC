package samplingnetwork

import (
	"fmt"
	"testing"

	"gonum.org/v1/gonum/mat"
)

func TestIntegratorBlock(t *testing.T) {
	gain := 1e2
	iBlock := IntegratorBlock(gain)
	if iBlock.System.InputSpaceOrder() != 1 {
		t.Error("Input Space Error")
	}

	if iBlock.System.OutputSpaceOrder() != 1 {
		t.Error("Output space error")
	}

	if iBlock.System.StateSpaceOrder() != 1 {
		t.Error("State space order error")
	}

	if len(iBlock.Control) > 1 {
		t.Error("Should only be one control")
	}

	iBlock.Control[0].SetState(0)
	if iBlock.Control[0].GetState() != uint(0) {
		t.Error("Get state is making an error")
	}

	a := iBlock.System.A.At(0, 0)
	b := iBlock.System.B.At(0, 0)
	c := iBlock.System.C.At(0, 0)

	if a != 0 || b != gain || c != 1 {
		t.Error("LinearSystem badly initialized")
	}
}

func TestSplitBlock(t *testing.T) {
	gain := 10.
	i1 := IntegratorBlock(gain)
	i2 := IntegratorBlock(gain)
	o1 := OscillatorBlock(gain, 1)
	o2 := OscillatorBlock(gain, 1)
	o3 := OscillatorBlock(gain, 1)

	sB1 := SplitBlock([]SamplingNetwork{i1, i2})
	sB2 := SplitBlock([]SamplingNetwork{o1, o2, o3})

	fmt.Println(mat.Formatted(sB2.Control[1].GetVector()))
	fmt.Printf("A = \n%v\nB = \n%v\nC =\n%v\n", mat.Formatted(sB1.System.A), mat.Formatted(sB1.System.B), mat.Formatted(sB1.System.C))
	fmt.Printf("A = \n%v\nB = \n%v\nC =\n%v\n", mat.Formatted(sB2.System.A), mat.Formatted(sB2.System.B), mat.Formatted(sB2.System.C))
}

func TestSplitBlockPanic(t *testing.T) {
	defer func() {
		if r := recover(); r == nil {
			t.Error("The code did not panic for different input dimensions")
		}
		fmt.Println("Code panic worked for split block.")
	}()
	gain := 1.
	i := IntegratorBlock(gain)
	o := OscillatorBlock(gain, gain)

	SplitBlock([]SamplingNetwork{i, o})
}

func TestMergeBlock(t *testing.T) {
	gain := 1.
	resonanceFrequency := 10.
	o1 := OscillatorBlock(gain, resonanceFrequency)
	o2 := OscillatorBlock(gain, resonanceFrequency)

	sys := MergeBlock([]SamplingNetwork{o1, o2})
	fmt.Printf("A = \n%v\nB = \n%v\nC =\n%v\n", mat.Formatted(sys.System.A), mat.Formatted(sys.System.B), mat.Formatted(sys.System.C))
}

func TestParallelBlock(t *testing.T) {
	gain := 1.
	resonanceFrequency := 10.
	o1 := OscillatorBlock(gain, resonanceFrequency)
	o2 := OscillatorBlock(gain, resonanceFrequency)

	sys := ParallelBlock([]SamplingNetwork{o1, o2})
	fmt.Printf("A = \n%v\nB = \n%v\nC =\n%v\n", mat.Formatted(sys.System.A), mat.Formatted(sys.System.B), mat.Formatted(sys.System.C))
}

func TestSeriesBlock(t *testing.T) {
	gain := 1.
	resonanceFrequency := 10.
	o1 := OscillatorBlock(gain, resonanceFrequency/1.)
	o2 := OscillatorBlock(gain*2, resonanceFrequency/10.)
	o3 := OscillatorBlock(gain*3, resonanceFrequency/1e2)
	o4 := OscillatorBlock(gain*4, resonanceFrequency/1e3)
	o5 := OscillatorBlock(gain*5, resonanceFrequency/1e4)

	sys := SeriesBlock([]SamplingNetwork{o1, o2, o3, o4, o5})
	fmt.Printf("A = \n%v\nB = \n%v\nC =\n%v\n", mat.Formatted(sys.System.A), mat.Formatted(sys.System.B), mat.Formatted(sys.System.C))
}

func TestIntegratorChain(t *testing.T) {
	gain := 6250.
	N := 5
	var integrators []SamplingNetwork
	for index := 0; index < N; index++ {
		integrators = append(integrators, IntegratorBlock(gain))
	}

	sys := SeriesBlock(integrators)
	fmt.Printf("A = \n%v\nB = \n%v\nC =\n%v\n", mat.Formatted(sys.System.A), mat.Formatted(sys.System.B), mat.Formatted(sys.System.C))
}
