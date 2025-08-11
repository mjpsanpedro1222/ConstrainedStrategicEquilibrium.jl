# # Symmetric - asymmetric CSE test case
#
# This example is using the asymmetric afrprogs code with two players with the
# same distribution (i.e. symmetric). The purpose is to test that the asymmetric
# code base can reproduce the symmetric result and we can also test this further
# by comparing the result to the BNE.
#
# ## Install dependencies
#
# Enter the Julia REPL (`julia`) and run:
#
# ```julia
# using Pkg
# Pkg.add("https://github.com/mjpsanpedro1222/ConstrainedStrategicEquilibrium.jl")
# ```

# ## Load required modules

using ConstrainedStrategicEquilibrium
using Plots
using Distributions
using NonlinearSolve

# ## Define some parameters
#
# Here we set the initial and final values of n.

inin = 2
maxn = 6

# ## Solve for 2 players
#
# Start from n=2 like the original sym.f code with the following initial guess:

xguess = [1.58285, 1.02813, -1.66483]

# Create the problem:

cse_prob2 = AsymmetricAfrprogsCSEProblem(
    inin=inin,
    maxn=maxn,
    np=2,
    mc=10000,
    solver_kwargs=(show_trace=Val(false), maxiters=2000, abstol=1e-6, reltol=1e-6),
    solver_initial_guess=xguess,
    distributions=[Beta(3, 3), Beta(3, 3)],
    knot_refinement_strategy=:even_spacing,
)

# Compute the solution to the CSE problem for 2 players:

solutions2 = compute_cse(cse_prob2)

# We also compute the BNE for 2 players for the same problem:

bne2x = Vector{Float64}(undef, 101)
bne2y = Vector{Float64}(undef, 101)
for m = 1:101
    ti = (m - 1.0) / 100.0
    if ti == 0
        true_bne = 0.0
    else
        true_bne = compute_bne(ti, cse_prob2.distributions[1], cse_prob2.np)
    end
    bne2x[m] = ti
    bne2y[m] = true_bne
end


# ## Solve for 4 players
#
# Define the initial guess and create the problem:
xguess = [-log(0.25), -log(0.25), log(0.25)]
cse_prob4 = AsymmetricAfrprogsCSEProblem(
    inin=inin,
    maxn=maxn,
    np=4,
    mc=10000,
    solver_kwargs=(show_trace=Val(false), maxiters=2000, abstol=1e-6, reltol=1e-6),
    solver_initial_guess=xguess,
    distributions=[Beta(3, 3), Beta(3, 3)],
    knot_refinement_strategy=:even_spacing,
)

# Compute the solution to the CSE problem for 4 players:

solutions4 = compute_cse(cse_prob4)

# We also compute the BNE for 4 players for the same problem:

bne4x = Vector{Float64}(undef, 101)
bne4y = Vector{Float64}(undef, 101)
for m = 1:101
    ti = (m - 1.0) / 100.0
    if ti == 0
        true_bne = 0.0
    else
        true_bne = compute_bne(ti, cse_prob4.distributions[1], cse_prob4.np)
    end
    bne4x[m] = ti
    bne4y[m] = true_bne
end

# ## Results
#
# First some sanity checking that both problems were solved successfully.

if length(solutions2) != length(solutions2)
    throw("solutions are not the same length")
end
if (!solutions2[end].success) || (!solutions4[end].success)
    throw("final solution not successful")
end

# Plot the two CSE solutions and BNE solutions on the same graph:

cseplot(
    solutions2[end],
    cse_label=Dict(:bidder1 => "CSE 1 (np=2)", :bidder2 => "CSE 2 (np=2)"),
    cse_colour=Dict(:bidder1 => :1, :bidder2 => :2),
    knot_label=Dict(:bidder1 => "", :bidder2 => ""),
    dpi=300,
)
cseplot!(
    solutions4[end],
    cse_label=Dict(:bidder1 => "CSE 1 (np=4)", :bidder2 => "CSE 2 (np=4)"),
    cse_colour=Dict(:bidder1 => :3, :bidder2 => :4),
    knot_label=Dict(:bidder1 => "", :bidder2 => ""),
    title="Comparing CSE against BNE for num players 2 and 4\nasym code but all bidders are Beta(3,3)",
    dpi=300,
)
plot!(bne2x, bne2y, label="BNE (np=2)", dpi=300, seriescolor=:5)
plot!(bne4x, bne4y, label="BNE (np=4)", dpi=300, seriescolor=:6)

## save the figure to file: #jl
savefig("afr-progs-asym-sym-compare.png")  #jl
