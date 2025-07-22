# # Symmetric CSE example - polynomial form
#
# This is a simple example of how to run the jae\_poly\_1 symmetric CSE case from Computer\_Code\_CSE.
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
using NonlinearSolve
using Distributions

# ## Create the symmetric CSE problem
#
# Create the default symmetric CSE problem with the following:

cse_prob = SymmetricJaePoly1CSEProblem()

# You can view the default options by running:

Base.doc(SymmetricJaePoly1CSEProblem)

# ## Compute the CSE
#
# Now compute the CSE for the problem we created
solutions = compute_cse(cse_prob)

# ## Postprocessing
#
# Plot the final successful solution comparing the computed CSE to the "analytical" Bayes-Nash Equilibrium.

for sol in Iterators.reverse(solutions)
    if sol.success
        plot(sol, dpi=400)
        savefig("jae-poly-1.png")
        break
    end
end

# View the plot showing the CSE and BNE: ![jae-poly-1.png](jae-poly-1.png)

# ## Create a different problem and solve it
#
# Now create a non-default problem, e.g. by changing the number of players, the distribution
# and the maximum value of n:

cse_prob = SymmetricJaePoly1CSEProblem(np=2, distribution=Beta(3, 3), maxn=14)

# Compute the CSE:

solutions = compute_cse(cse_prob)

# Finally, plot the last successful solution:

for sol in Iterators.reverse(solutions)
    if sol.success
        plot(sol, dpi=400)
        savefig("jae-poly-1-beta.png")
        break
    end
end

# View the plot showing the CSE and BNE: ![jae-poly-1-beta.png](jae-poly-1-beta.png)
