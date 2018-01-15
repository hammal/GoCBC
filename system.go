package adc

// System struct contains all relevant system parameters for the simulator
type System struct {
	// Time period
	Ts float64
	// starting time
	StartTime float64
	// ending time
	EndTime float64
	// Number of samples
	N int64
	// Number of Inputs
	NumberOfInputs int
	// Number of observations
	NumberOfObservations int
}
