# ConstrainedStrategicEquilibrium.jl

ConstrainedStrategicEquilibrium.jl is a Julia package for solving Constrained
Strategic Equilibrium (CSE) problems in game theory and economics.

It has been implemented based on the work by Armantier et al [armantier2008cse](@cite)
and the Fortran code released alongside that publication.

## Features

- **Multiple Problem Types**: Symmetric and asymmetric CSE problems
- **Flexible Algorithms**: Support for NonlinearSolve.jl solvers
- **Built-in Visualization**: Plotting recipes with RecipesBase
- **Data Export**: Results as DataFrames for easy analysis

## Problem Types Supported

- [`SymmetricAfrprogsCSEProblem`](@ref) - Symmetric piecewise linear CSE form
- [`SymmetricJaePoly1CSEProblem`](@ref) - Symmetric polynomial CSE form
- [`AsymmetricAfrprogsCSEProblem`](@ref) - Asymmetric piecewise linear CSE form

## Quick start

Install the CSE package:

```
julia -e 'using Pkg; Pkg.add(url="https://github.com/mjpsanpedro1222/ConstrainedStrategicEquilibrium.jl")'
```

Then in the julia repl (enter by running `julia`):

```julia
using ConstrainedStrategicEquilbrium

# create a CSE problem with default options
prob = SymmetricJaePoly1CSEProblem()

# solve the CSE problem
solutions = compute_cse(prob)

# print the vector of solutions
println(solutions)

# plot the last solution
using Plots
cseplot(solutions[end]; dpi=300)
savefig("jae-poly-1-result.png")

# store the computed CSE/BNE to csv file
using CSV
CSV.write("cse_result.csv", solutions[end].cse)
```
