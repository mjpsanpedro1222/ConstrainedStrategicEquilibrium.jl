# Contrained Strategic Equilibrium Julia package

ConstrainedStrategicEquilibrium.jl is a Julia package for solving Constrained
Strategic Equilibrium (CSE) problems in game theory and economics.

It has been implemented based on the work by Armantier et al and the Fortran
code released alongside their paper [https://doi.org/10.1002/jae.1040](https://doi.org/10.1002/jae.1040).

Additional documentation can be found [here](https://mjpsanpedro1222.github.io/ConstrainedStrategicEquilibrium.jl/).

## Install dependencies

```
julia -e 'using Pkg; Pkg.add("https://github.com/mjpsanpedro1222/ConstrainedStrategicEquilibrium.jl")'
```

## Run a simple case

```julia
using ConstrainedStrategicEquilibrium

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

## Run an example

Enter the Julia REPL (for example, by running `julia`).

Run the example symmetric calculation (from afr-progs) using:

```julia
julia> include("examples/afr-progs-sym.jl")
```

You could then investigate the solutions for different `n` values, e.g.

```julia
# display the vector of solutions
julia> display(solutions)

# select the last solution
julia> sol = solutions[end]

# plot it
julia> plot(sol)

# write the CSE/BNE to CSV file
julia> using CSV
julia> CSV.write("cse_result.csv", sol.cse)

# this can be read back in with
julia> using DataFrames
julia> df = CSV.read("cse_result.csv", DataFrame)
```

Many of the other examples can also be run this way, for example:

```julia
julia> include("examples/jae-poly-1.jl")
```

## Build the docs

You can build the documentation locally, which can be useful while developing the docs further.

First, you need to instantiate the docs project:

```
julia --project=docs -e "using Pkg; Pkg.instantiate()"
```

Then build the docs:

```
julia --project=docs docs/make.jl
```

Open *docs/build/index.html* in your browser to view the local version of the docs.

## Enable debug logging

We use the Logging package to print info. To enable debug messages you could
export this environment variable before running julia:

```
export JULIA_DEBUG=ConstrainedStrategicEquilibrium
```
