# State Space Model Implementation in Go
Author: Hampus Malmberg

This is my state space implementation of Go.

## Overview
This package is fundamentally intended to hold a state space model specified by a
- State transition matrix A
- Input vector B
- Observation Vector C
- and D vector

Additionally, it should provide the functions
- Step(output, input) which computes the next step as $X = AX + BU$
- Observation() ....
