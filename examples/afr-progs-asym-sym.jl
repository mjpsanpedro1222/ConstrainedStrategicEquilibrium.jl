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
    maxn=12,
    np=4,
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
# Compute the BNE

bnex = Vector{Float64}(undef, 101)
bney = Vector{Float64}(undef, 101)
for m = 1:101
    ti = (m - 1.0) / 100.0
    if ti == 0
        true_bne = 0.0
    else
        true_bne = compute_bne(ti, cse_prob.distributions[1], cse_prob.np)
    end
    bnex[m] = ti
    bney[m] = true_bne
    println("BNE($(ti)) = $(true_bne)")
end

#
#
# Loop over the solutions and plot them all:
for sol in solutions
    if sol.success
        println("Plotting n=$(sol.n)")

        # compute MSE against BNE
        mse1_sum = 0.0
        mse2_sum = 0.0
        for m = 1:101
            if m == 1
                cse1 = cse2 = 0
            else
                cse1 = sol.cse[!, "CSE(x) 1"][m]
                cse2 = sol.cse[!, "CSE(x) 2"][m]
            end
            #println("CSEs: $(cse1); $(cse2)")
            mse1_sum += (sol.cse[!, "CSE(x) 1"][m] - bney[m])^2
            mse2_sum += (sol.cse[!, "CSE(x) 2"][m] - bney[m])^2
        end
        mse1 = mse1_sum / 101
        mse2 = mse2_sum / 101
        println("MSE bidder 1: $(mse1)")
        println("MSE bidder 2: $(mse2)")

        cseplot(sol, dpi=300)
        plot!(bnex, bney, label="BNE")

        # We can also save the figure to a file (comment this line out if you don't want to)
        savefig("afr-progs-asym-n$(sol.n).png")
    end
end
