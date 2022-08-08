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

## References
D. Xiu and G. E. Karniadakis, “Modeling uncertainty in flow simulations via Generalized Polynomial Chaos,” Journal of Computational Physics, vol. 187, no. 1, pp. 137–167, 2003.

D. Xiu and G. E. Karniadakis, “The wiener--askey polynomial chaos for stochastic differential equations,” SIAM Journal on Scientific Computing, vol. 24, no. 2, pp. 619–644, 2002.

Y.-zhong Luo and Z. Yang, “A review of uncertainty propagation in Orbital Mechanics,” Progress in Aerospace Sciences, vol. 89, pp. 23–39, 2017.
