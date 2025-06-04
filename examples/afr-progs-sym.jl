# # Symmetric CSE example - piecewise linear form
#
# This is a simple example of how to run the symmetric CSE case from afr_progs.
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
# Select the final solution from the list
sol = solutions[end]

# Plot the final solution comparing the computed CSE to the "analytical" Bayes-Nash Equilibrium and save the figure to a file.
plot(sol, dpi=300)
savefig("afr_progs_sym.png")

# View the plot showing the CSE and BNE: ![afr_progs_sym.png](afr_progs_sym.png)

# ## Create a different problem and solve it
#
# Now create a non-default problem, e.g. by lowering the maximum value of n and changing the distribution parameters:

cse_prob = SymmetricAfrprogsCSEProblem(maxn=12, distribution=Beta(2.5, 3.5))

# Compute the CSE:

solutions = compute_cse(cse_prob)

# Finally, plot the last solution:

plot(solutions[end], dpi=400)
savefig("afr_progs_sym_n12.png")

# View the plot showing the CSE and BNE: ![afr_progs_sym_n12.png](afr_progs_sym_n12.png)
