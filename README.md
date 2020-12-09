# BioSANS - Symbolic and Numeric Software for Systems Biology
*BioSANS* is a free software for systems biology which is currently being developed in Academia Sinica Institute of Chemistry. The goal of this development is to make systems biology available to non-domain expersts in the easiest possible way. The software support SBML and passed majority of the semantic and stochastic test cases. We also tested the software's capability to perform symbolic computation (species analytical concentration) using 40 symbolic test cases. This test cases where used as well to test parameter estimation capability with additional 5 cases for ill conditioned parameters order of magnitude. Symbolic computation for linear noise approximation was tested against 20 test cases that can go to steady state. 

The following summarized the symbolic and numeric features currently supported;

### Symbolic computation

1. Species analytical expression - works for most linear differential equation and few non-linear ordinary differential equations
2. LNA covariace matrix - works for most linear differential equation and few non-linear ODE
3. Steady state concentration - generally works for most problems especially linear ODE
4. Network localization - topology based sentitivity matrix

### Numeric computation

1. Linear noise approximation
2. Parameter estimation
3. Network localization
4. Deterministic analysis (ODE integration)
  * odeint - pythpn LSODA library
  * runge-kutta4 (tau-adaptive and fixed interval version)
  * Euler (two differnt types of tau-adaptive version)
5. Stochastic modeling
  * Chemical Langevine equation (tau-adaptive and fixed interval version)
  * Tau-leaping algorithm (3 different versions/implementation of Yang Cao's algorithm)
  * Gillespie direct method
  
### Post-processing functions

1. Plotting (trajectory, density, etc.)
2. Calculation of correlation, covariance, fano-factor, etc.
3. etc.
