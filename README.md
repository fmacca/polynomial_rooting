# polynomials-simulations
Simulations on roots of perturbed polynomials

This repository contains:
- Simulation scripts, with format `DegreeN_NoiseXxxxx_Yyyyyyy`, where:
    - `N`: Degree of the polynomial
    - `Xxxxx`: Circular, FullMatrix
    - `Yyyyyyy`: a description of the parameter that is simulated
    This are the simulations described in Section 5 Numerical Analysis of the paper.
- Logistic regression scripts, computing thresholds of validity for the simulated parameter.
    They load the data produced by the simulation scripts.
- Other simulations, reproducing the other figures in the paper.

## Requirements
- Matlab
- Statistics and Machine Learning Toolbox
