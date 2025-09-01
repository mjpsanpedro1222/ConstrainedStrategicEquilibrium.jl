# # Hyperparameter optimisation for initial conditions
#
# This example demonstrates the use of HyperTuning to find good initial conditions
# for an asymmetric auction problem with 2 players. Note this example uses the 
# asymmetric code base but both players have the same distribution so is effectively
# symmetric.
#
# ## Install dependencies
#
# Enter the Julia REPL (`julia`) and run:
#
# ```julia
# using Pkg
# Pkg.add("ConstrainedStrategicEquilibrium")
# Pkg.add("HyperTuning")
# ```
#
# ## Load required modules

using ConstrainedStrategicEquilibrium
using Plots
using Distributions
using NonlinearSolve
using HyperTuning
using Logging

# ## Define the objective function
#
# The objective function solves the CSE and computes the Mean Squared Error (MSE) against the Bayesian
# Nash Equilibrium (BNE). It takes trial values for a, b, c and uses them as initial guesses.
#
# We are able to compute the BNE (and therefore MSE) because both players have the same distribution. If you
# were to run this for a true asymmetric case then you would need to choose another metric to return
# from the `objective` function, such as `sol.c_2` or the norm of the residual vector from the solver.
#
# We create the CSE problem only for the initial value of `n` (2 in this case) and set the tolerances loosely
# as we only need to roughly converge on a good initial guess.

function objective(trial)
    ## unpack the initial guess and store it in a vector suitable for passing to the problem
    @unpack a, b, c = trial
    xguess = [a, b, c]

    ## create the problem setting inin==maxn
    cse_prob = AsymmetricAfrprogsCSEProblem(
        inin=2,
        maxn=2,
        np=2,
        mc=10000,
        solver_kwargs=(show_trace=Val(false), maxiters=2000, abstol=1e-6, reltol=1e-6),
        solver_initial_guess=xguess,
        distributions=[Beta(3, 3), Beta(3, 3)],
        knot_refinement_strategy=:even_spacing,
    )

    ## solve the CSE
    solutions = compute_cse(cse_prob)

    ## now we compute the BNE for the same problem
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
    end

    ## compute the MSE for the CSE against the BNE and return that
    mse = Inf
    sol = solutions[1]
    if sol.success
        mse1_sum = 0.0
        mse2_sum = 0.0
        for m = 1:101
            if m == 1
                cse1 = cse2 = 0
            else
                cse1 = sol.cse[!, "CSE(x) 1"][m]
                cse2 = sol.cse[!, "CSE(x) 2"][m]
            end
            mse1_sum += (sol.cse[!, "CSE(x) 1"][m] - bney[m])^2
            mse2_sum += (sol.cse[!, "CSE(x) 2"][m] - bney[m])^2
        end
        mse1 = mse1_sum / 101
        mse2 = mse2_sum / 101
        mse = mse1 + mse2
    end

    return mse
end

# ## Run the optimisation
#
# Create the HyperTuning `Scenario` and define some ranges for the initial guess
# parameters (these could be made wider if needed).

scenario = Scenario(
    a=(0.1 .. 3.0),
    b=(0.1 .. 3.0),
    c=(-3.0 .. -0.1),
    max_trials=400,
)

# Run the optimisation, noting that we set the log level to be really high so that we
# don't get inundated with messages from the hundreds of trials being run, and print
# the best initial guess that was found.

logger = ConsoleLogger(stderr, Logging.AboveMaxLevel)
best_guess = with_logger(logger) do
    sc = HyperTuning.optimize(objective, scenario)
    println(sc)  #jl
    best_guess = [sc.best_trial.values[:a], sc.best_trial.values[:b], sc.best_trial.values[:c]]

    return best_guess
end
println("Best guess initial condition: $best_guess")

# ## Use the best guess
#
# Create a full problem, in this case solving for `n=2..16`, starting from the best initial
# guess that we found above.

cse_prob = AsymmetricAfrprogsCSEProblem(
    inin=2,
    maxn=16,
    np=2,
    mc=10000,
    solver_kwargs=(show_trace=Val(false), maxiters=2000, abstol=1e-6),
    solver_initial_guess=best_guess,
    distributions=[Beta(3, 3), Beta(3, 3)],
    knot_refinement_strategy=:even_spacing,
)

# Solve the problem.

solutions = compute_cse(cse_prob)

# Compute the BNE too

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
end

# Plot the last CSE solution and the BNE

cseplot(solutions[end])
plot!(bnex, bney, label="BNE")

## save the figure to file: #jl
savefig("afr-progs-asym-sym-compare.png")  #jl
