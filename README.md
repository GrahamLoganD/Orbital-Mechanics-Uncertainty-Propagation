# Orbital Mechanics Uncertainty Propagation

These programs propagate the uncertainty in a satellite's state vector (i.e. its position and velocity) over time as it orbits the Earth. This has application in space situational awareness and trajectory design/optimization. The programs are written in Python and Matlab, using Monte Carlo and Polynomial Chaos Expansion methods. An implementation of Polynomial Chaos Expansion in the UQLab framework for uncertainty quantification is included as well.

**Warning:** These programs can be computationally intensive and take over five minutes to run, depending on the settings.

### Matlab
- [Monte Carlo Simulation](https://github.com/GrahamLoganD/Orbital-Mechanics-Uncertainty-Propagation/blob/master/MC_Simulation.m)
- [Polynomial Chaos Expansion Simulation](https://github.com/GrahamLoganD/Orbital-Mechanics-Uncertainty-Propagation/blob/master/PCE_Simulation.m)
- [Polynomial Chaos Expansion Simulation using UQLab](https://github.com/GrahamLoganD/Orbital-Mechanics-Uncertainty-Propagation/blob/master/UQLab%20PCE%20Simulation/UQLab_PCE_Simulation.m)

### Python
- [Monte Carlo Simulation](https://github.com/GrahamLoganD/Orbital-Mechanics-Uncertainty-Propagation/blob/master/MC_Simulation.py)
- [Polynomial Chaos Expansion Simulation](https://github.com/GrahamLoganD/Orbital-Mechanics-Uncertainty-Propagation/blob/master/PCE_Simulation.py)

## Dependencies
- [UQLab](https://www.uqlab.com/)
- numpy
- scipy
- pandas
- seaborn
- matplotlib
