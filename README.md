# Contrained Strategic Equilibrium Julia package

Prototype of how we could have a Julia package... would be best to move this to its own repo

## Install dependencies

From this directory run:

```
julia --project
```

then run:

```julia
julia> using Pkg
julia> Pkg.instantiate()
```

## Run symmetric example

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

You could also run the symmetric jae_poly_1 example:

```julia
julia> include("examples/jae-poly-1.jl")
```

and so on...

## Build the docs

Need to instantiate the docs project first:

```
julia --project=docs -e "using Pkg; Pkg.instantiate()"
```

Then build the docs:

```
julia --project=docs docs/make.jl
```

Open *docs/build/index.html*.

## Enable debug logging

We use the Logging package to print info. To enable debug messages you could
export this environment variable before running julia:

```
export JULIA_DEBUG=ConstrainedStrategicEquilibrium
```

## More to come...
