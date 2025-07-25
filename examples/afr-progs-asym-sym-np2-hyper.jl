# this example using HyperTuning to try to find the best initial conditions

# ## Load required modules

using ConstrainedStrategicEquilibrium
using Plots
using Distributions
using NonlinearSolve
using HyperTuning


function objective(trial)
    @unpack a, b, c = trial
    xguess = [a, b, c]

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

    validate_cse_problem(cse_prob)

    solutions = compute_cse(cse_prob)

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

    mse = Inf
    sol = solutions[1]
    if sol.success
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
            mse1_sum += (sol.cse[!, "CSE(x) 1"][m] - bney[m])^2
            mse2_sum += (sol.cse[!, "CSE(x) 2"][m] - bney[m])^2
        end
        mse1 = mse1_sum / 101
        mse2 = mse2_sum / 101
        mse = mse1 + mse2
    end

    return mse
end



function main()
    scenario = Scenario(
        a=(0.1 .. 3.0),
        b=(0.1 .. 3.0),
        c=(-3.0 .. -0.1),
        max_trials=400,
    )

    sc = HyperTuning.optimize(objective, scenario)
    println("RESULT:")
    println(sc)
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
