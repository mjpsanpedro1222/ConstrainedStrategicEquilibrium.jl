# # Symmetric CSE example
#
# This is a simple example of how to run the jae_poly_1 symmetric CSE case from Compute_Code_CSE.
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
# Select the final solution from the list
sol = solutions[end]

# Plot the final solution comparing the computed CSE to the "analytical" Bayes-Nash Equilibrium.
plot(sol, dpi=400)
savefig("jae_poly_1.png")

# ## Create a different problem and solve it
#
# Now create a non-default problem, e.g. by increasing the value of n and the number of players:

cse_prob = SymmetricJaePoly1CSEProblem(np = 6, maxn = 9)

# Compute the CSE:

solutions = compute_cse(cse_prob)

# Finally, plot the last solution:

plot(solutions[end], dpi = 400)
savefig("jae_poly_1_n9.png")
