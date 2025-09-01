# # Symmetric CSE example - polynomial form
#
# This is a simple example of how to run the jae\_poly\_1 symmetric CSE case from Computer\_Code\_CSE
# coming from the work by Armantier et al. [armantier2008cse](@cite).
#
# ## Install CSE package
#
# Enter the Julia REPL (`julia`) and run:
#
# ```julia
# using Pkg
# Pkg.add("ConstrainedStrategicEquilibrium")
# ```

# ## Load required modules

using ConstrainedStrategicEquilibrium
using Plots
using NonlinearSolve
using Distributions

# ## Create the symmetric CSE problem
#
# Create the CSE problem with two players and a maximum n value of 12:

cse_prob = SymmetricJaePoly1CSEProblem(np=2, maxn=12)

# ## Compute the CSE
#
# Now compute the CSE for the problem we created
solutions = compute_cse(cse_prob)

# ## Postprocessing
#
# Plot the final successful solution comparing the computed CSE to the "analytical" Bayes-Nash Equilibrium.

for sol in Iterators.reverse(solutions)
    if sol.success
        cseplot(sol, dpi=400)
        savefig("jae-poly-1.png")
        break
    end
end

# View the plot showing the CSE and BNE: ![jae-poly-1.png](jae-poly-1.png)

# ## References
#
# * [armantier2008cse](@cite) Armantier et al. Journal of Applied Econometrics, 23 (2008)
