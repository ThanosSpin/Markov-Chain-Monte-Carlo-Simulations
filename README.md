# Markov-Chain-Monte-Carlo-Simulations
Markov Chain Monte Carlo (MCMC) Simulations for Volume Estimation of a d-dimensional Sphere


The goal of this repository is the use of MCMC Simulations for volume estimation of a d-dimensional Sphere.

To begin with, I use MCMC in order to estimate the integral of a circle in two dimensions with radius R = 1. 
Then, I repeat the process for more dimensions (d = 1,.., d = 200) and I compare the results with MC simulations. 

The volume can be estimated by running a Markov Chain inside a high - dimensional sphere (Sd-1). The Energy of Metropolis Hastings is set infinity when x not in Sd and 1 when x in Sd. I run the chain in Sd-1 and then I add a dimension independently in [-1,1] and the test is if the simulated values x1^2 + x2^2 + ... + xd^2 are <= 1 (so the chain remains in Sd sphere) or > 1 so the chain after the added dimension has moved to cylinder, d - dimensions, for R = 1.

Knowing that P = |Sd|/|Cyld|=|Sd|/2*|Sd-1| the volume of |Sd| can be estimated: |Sd| = 2*P*|Sd-1| (|Sd|: Volume of Sphere d-dimensions, |Cyld|: Volume of |Cylinder d-dimensions|) 
