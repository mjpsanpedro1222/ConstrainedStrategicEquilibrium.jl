# # Symmetric CSE example - piecewise linear form
#
# This is a simple example of how to run the symmetric CSE case from the supporting material
# to the paper by Armantier et al. [armantier2008cse](@cite).
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

# ## Create the symmetric CSE problem
#
# Create the default symmetric CSE problem with the following:

cse_prob = SymmetricAfrprogsCSEProblem()

# ## Compute the CSE
#
# Now compute the CSE for the problem we created
solutions = compute_cse(cse_prob)

# ## Postprocessing
#
# Plot the final successful solution, comparing the computed CSE to the "analytical"
# Bayes-Nash Equilibrium and save the figure to a file.

for sol in Iterators.reverse(solutions)
    if sol.success
        cseplot(sol; dpi=300)
        savefig("afr-progs-sym.png")
        break
    end
end

# View the plot showing the CSE and BNE: ![afr-progs-sym.png](afr-progs-sym.png)

# ## Create a different problem and solve it
#
# Now create a non-default problem, e.g. by lowering the maximum value of n and changing the distribution parameters:

cse_prob = SymmetricAfrprogsCSEProblem(maxn=12, distribution=Beta(2.5, 3.5))

# Compute the CSE:

solutions = compute_cse(cse_prob)

# Finally, plot the last successful solution:

for sol in Iterators.reverse(solutions)
    if sol.success
        cseplot(sol; dpi=400)
        savefig("afr-progs-sym-beta2535.png")
        break
    end
end

# View the plot showing the CSE and BNE: ![afr-progs-sym-beta2535.png](afr-progs-sym-beta2535.png)

# ## References
#
# * [armantier2008cse](@cite) Armantier et al. Journal of Applied Econometrics, 23 (2008)
