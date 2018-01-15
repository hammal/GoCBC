package control

// Control interface holds the...
type Control interface {
	// Executes the control simulation
	Simulate()
	// Get control sequence at index
	Get(index int) []bool
}
