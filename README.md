# Analog to digital converter simulator in Go
This is a Go library for the analog to digital (ADC) concepts concerning my PhD.

## Code Structure
The main components of this implementation is the:
#### ADC
which is the main interface and in-turn orchestrate the overall actions such as:
  - Simulate(), which initiates an simulation on the current instance.
  - Reconstruct(), which uses the control information to make an reconstruction of signal/signals that where feed into the network.
  - Load(), Load a control sequence into this instance.
  - Save(), Save the current instance for later use.
  - GetTimeStamps(), return the absolute time stamps for the internal indexes.

Additionally, this instance keeps the necessary types to describe the relevant system
information such as time constants, framework configurations and absolute quantities.


#### Control
This instance holds the current control decisions, the basis for reconstruction, and is one of
the main components for simulating the ADC framework.

#### Simulator
The main component for simulating an ADC framework requires an Control instance and
one or several input signals for simulation.

#### Reconstruction
This instance can with the help of an Control instance perform a full reconstruction.

## File Structure
This program has one main module that is used for interacting. These are
situated in the adc.go file Interface ADC. Furthermore the standard ADC
are created using one of the functions named New_... located in the same file. Additionally, there are some helper modules:
- [control](control/README.md), implements the control object
- [ode](ode/README.md), a helper module for doing standard ode solving.
- [reconstruct](reconstruct/README.md), implements the reconstruction framework.
- [signal](signal/README.md), implements the different signal types
- [simulate](simulate/README.md), implements the simulator which is used to simulate the ADC network.
- [ssm](ssm/README.md), implements different state space model formulations


## Notes
- ~~implement control pre computations~~
- TODO: implement adc general Interface
- ~~implement reconstruction steady state computations~~
- IDEA: Check out [Cobra](https://github.com/spf13/cobra) for command line flags
- IDEA: Check out [Viper](https://github.com/spf13/viper) for configuration files
- TODO: Post-filteing implementation
- TODO: Oscillator Control
- TODO: Eigenvalue decomposition for reconstruction 
