# # Symmetric CSE from afr-progs using the asymmetric code
#
# This example is using the asymmetric afrprogs code with two players with the
# same distribution (so symmetric).
#
# ## Install dependencies
#
# Make sure the required packages are installed by running the following from
# the root directory of this repo:

# ```julia
# using Pkg
# Pkg.instantiate()
# ```

# ## Load required modules

using ConstrainedStrategicEquilibrium
using Plots
using Distributions
using NonlinearSolve
using MINPACK

# ## Create the CSE problem
#
# Start from n=2 like the original sym.f code with the following initial guess that seems
# to work OK initially:

nval = 2
xguess = Vector{Float64}(undef, 2 * nval - 1)
xguess[1] = -log(0.25)
xguess[2] = -log(0.25)
xguess[3] = log(0.25)

#
# Create an asymmetric CSE problem with two players and both distributions the same:
cse_prob = AsymmetricAfrprogsCSEProblem(
    inin=nval,
    maxn=7,
    np=2,
    mc=10000,
    solver_kwargs=(show_trace=Val(true), maxiters=2000, abstol=1e-6, reltol=1e-6),
    solver_initial_guess=xguess,
    distributions=[Beta(3, 3), Beta(3, 3)],
    knot_refinement_strategy=:even_spacing,
)

# Now validate the problem we just created (this will throw an error if there is a mistake):
validate_cse_problem(cse_prob)

# ## Compute the CSE
#
# Now compute the CSE for the problem we created
solutions = compute_cse(cse_prob)

# ## Postprocessing
#
# Loop over the solutions and plot them all:
for sol in solutions
    if sol.success
        plot(sol, dpi=300)

        # We can also save the figure to a file (comment this line out if you don't want to)
        savefig("afr-progs-asym-n$(sol.n).png")
    end
end
