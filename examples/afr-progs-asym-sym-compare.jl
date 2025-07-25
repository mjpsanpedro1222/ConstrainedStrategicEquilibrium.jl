using Base: run_package_callbacks
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


function run_problem(cse_prob)
    # Now validate the problem we just created (this will throw an error if there is a mistake):
    validate_cse_problem(cse_prob)

    # ## Compute the CSE
    solutions = compute_cse(cse_prob)

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
    end

    return (solutions=solutions, bnex=bnex, bney=bney)
end


function main()
    inin = 2
    maxn = 6

    # ## Create the CSE problem for 2 players
    #
    # Start from n=2 like the original sym.f code with the following initial guess that seems
    # to work OK initially:

    xguess = [1.76034, 1.26109, -2.92289]
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

    # solve for 2 players
    result2 = run_problem(cse_prob2)

    ## now 4 players
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
    result4 = run_problem(cse_prob4)

    #
    # plot the last solution of each
    if length(result2.solutions) != length(result4.solutions)
        throw("solutions are not the same length")
    end
    if (!result2.solutions[end].success) || (!result4.solutions[end].success)
        throw("final solution not successful")
    end

    cseplot(
        result2.solutions[end],
        cse_label=Dict(:bidder1 => "CSE 1 (np=2)", :bidder2 => "CSE 2 (np=2)"),
        cse_colour=Dict(:bidder1 => :1, :bidder2 => :2),
        knot_label=Dict(:bidder1 => "", :bidder2 => ""),
        dpi=300,
    )

    cseplot!(
        result4.solutions[end],
        cse_label=Dict(:bidder1 => "CSE 1 (np=4)", :bidder2 => "CSE 2 (np=4)"),
        cse_colour=Dict(:bidder1 => :3, :bidder2 => :4),
        knot_label=Dict(:bidder1 => "", :bidder2 => ""),
        title="Comparing CSE against BNE for num players 2 and 4\nasym code but all bidders are Beta(3,3)",
        dpi=300,
    )

    plot!(result2.bnex, result2.bney, label="BNE (np=2)", dpi=300, seriescolor=:5)
    plot!(result4.bnex, result4.bney, label="BNE (np=4)", dpi=300, seriescolor=:6)

    savefig("afr-progs-asym-sym-compare.png")
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
