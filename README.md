# Strategic Policymaking for Implementing RPS: A Tri-level Optimization Approach

This code was developed under Julia v1.2.0/JuMP0.20 by Jip Kim in 2019.
The following packages must be installed:

  - JuMP
  - Distributions
  - LinearAlgebra
 
To run the code, execute CCG.jl, or include() it from a Julia prompt.

"ISONE8busTN" testsystem is given as follows in data folder:
  - Node.csv: Transmission network data
  - Line.csv: Transmission line data
  - Generator.csv: Generation data
  - We also specify additional input data in the beginning of MPini.jl/MP.jl/SP.jl
