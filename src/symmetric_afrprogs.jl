
# TODO: add stopping criteria tolerances; output files (data, results)
"""
$(TYPEDEF)

The symmetric CSE problem from afr-progs.

Parameters can be passed in as keyword arguments or can be omitted to accept the default values.

$(TYPEDFIELDS)

# Examples
```jldoctest
julia> prob = SymmetricAfrprogsCSEProblem()
SymmetricAfrprogsCSEProblem(np=4, mc=10000, n=2..16, Distributions.Beta{Float64}(α=3.0, β=3.0))

julia> prob = SymmetricAfrprogsCSEProblem(mc = 1000, maxn = 12, distribution = Beta(3, 4))
SymmetricAfrprogsCSEProblem(np=4, mc=1000, n=2..12, Distributions.Beta{Float64}(α=3.0, β=4.0))
```
"""
@kwdef struct SymmetricAfrprogsCSEProblem <: SymmetricCSEProblem
    "Random number generator to use during data generation (default rng is seeded with 642867)"
    rng::AbstractRNG = Random.seed!(642867)
    "Number of Monte Carlo steps (default is 10000)"
    mc::Int = 10000
    "Number of players (default is 4)"
    np::Int = 4
    "Distribution to use (must be `Beta` currently; default is `Beta(3, 3)`)"
    distribution::UnivariateDistribution = Beta(3, 3)
    "Initial value for n (default is 2 and must be left as this currently)"
    inin::Int = 2
    "Maximum value for n (default is 16)"
    maxn::Int = 16
    "Write txt and csv files with solution info"
    legacy_output::Bool = false
end


"""
$(TYPEDSIGNATURES)

Check that `cse_problem` is a valid definition for a CSE problem.

This function will throw an error if the problem is not valid, otherwise
will return silently.

# Examples
```jldoctest
julia> prob = SymmetricAfrprogsCSEProblem();
julia> validate_cse_problem(prob)

julia> prob = SymmetricAfrprogsCSEProblem(distribution = Normal());
julia> validate_cse_problem(prob)
ERROR: "Only Beta distributions are supported currently"
[...]

julia> prob = SymmetricAfrprogsCSEProblem(inin = 4);
julia> validate_cse_problem(prob)
ERROR: "Initial value of n must be 2 currently"
[...]

julia> prob = SymmetricAfrprogsCSEProblem(maxn = 1);
julia> validate_cse_problem(prob)
ERROR: "Initial value of n cannot be bigger than maximum value of n"
[...]
```
"""
function validate_cse_problem(cse_problem::SymmetricAfrprogsCSEProblem)
    if !(cse_problem.distribution isa Distributions.Beta)
        throw("Only Beta distributions are supported currently")
    end

    if cse_problem.np < 2
        throw("Not enough players")
    end

    if cse_problem.inin != 2
        throw("Initial value of n must be 2 currently")
    end

    if cse_problem.inin > cse_problem.maxn
        throw("Initial value of n cannot be bigger than maximum value of n")
    end
end


"""
$(TYPEDSIGNATURES)

Objective function for the symmetric afrprogs case.
"""
function objective_function_symmetric_afrprogs(fvec, x, p::SymmetricFunctionParams)
    # note: important to use similar here in case using autodiff they could be of type dual from ForwardDiff
    # TODO: preallocate for performance??
    da = similar(x)
    yknot = similar(x, length(x) + 1)
    alph = similar(da)
    bet = similar(da)

    # beta distribution
    # TODO: need to generalise for different distributions
    betadist = p.dist
    betadistparams = params(betadist)

    # set up the value of the constrained strategy parameters such that the
    # strategy is continuous
    # alph= constant, bet=slope, da(l)= derivative wrt l parameter

    # TODO: need to generalise for different distributions?
    const_val = beta(betadistparams...)

    da .= 0.0
    yknot[1] = 0.0
    for l = 1:p.n
        yknot[l+1] = yknot[l] + exp(x[l])
        alph[l] = yknot[l]
        bet[l] = (yknot[l+1] - yknot[l]) / (p.knot[l+1] - p.knot[l])
    end

    for m = 1:p.mc
        ti = p.u[m, 1]
        check = true
        l = 0
        while check
            l += 1
            if (ti >= p.knot[l]) && (ti <= p.knot[l+1])
                bi = alph[l] + bet[l] * (ti - p.knot[l])
                dbdt = bet[l]
                dbdp = (ti - p.knot[l]) / (p.knot[l+1] - p.knot[l])
                check = false

                # TODO: need to generalise for different distributions?
                cumu1 = cdf(betadist, ti)
                dcumu1 = (ti^(betadistparams[1] - 1)) * ((1 - ti)^(betadistparams[2] - 1)) / const_val
                cumu = cumu1^(p.np - 1)
                dcumu = (p.np - 1.0) * dcumu1 * cumu1^(p.np - 2)

                da[l] += dbdp * ((ti - bi) * dcumu / dbdt - cumu)
            end
        end
    end

    for l = 1:p.n
        fvec[l] = da[l] * yknot[l+1] / p.mc
    end

    # if cvrg flag is set, output results
    if p.cvrg
        # TODO: filenames should be an option
        # TODO: output everything back as a solution instead of writing to file here?
        if p.legacy_output
            fout = open("sym-result-n-$(p.n).txt", "w")
            fcsv = open("sym-bids-private-values-n-$(p.n).csv", "w")
            write(fout, "  t    CSE(t)    BNE\n")
            write(fcsv, "t, CSE(t), BNE\n")
        end

        mse_sum = 0.0
        count_valid = 0
        for m = 1:101
            ti = (m - 1.0) / 100.0
            bi = missing
            true_bne = compute_bne(ti, betadist, p.np)

            # If the first element's result is NaN, replace it with 0.0
            if m == 1 && isnan(true_bne)
                true_bne = 0.0
            end

            check = true
            l = 0
            while check
                l += 1
                if (ti >= p.knot[l] && ti <= p.knot[l+1])
                    bi = alph[l] + bet[l] * (ti - p.knot[l])
                    check = false
                end
            end

            # Compute squared error if bi is not missing
            if !ismissing(bi)
                mse_sum += (bi - true_bne)^2
                count_valid += 1
            end

            if p.legacy_output
                write(fout, "$ti    $bi    $true_bne\n")
                write(fcsv, "$ti, $bi, $true_bne\n")
            end
            push!(p.solution.cse, (ti, bi, true_bne))
        end

        if p.legacy_output
            close(fcsv)
        end

        # Compute MSE
        p.solution.mse = count_valid > 0 ? mse_sum / count_valid : NaN
        if p.legacy_output
            write(fout, "\nMean Squared Error: $(p.solution.mse)\n")

            write(fout, " \n l    knot(l-1)    knot(l)    CSE[knot(l)]\n")
        end
        for l = 1:p.n-1
            if p.legacy_output
                write(fout, "$(l*1.0) $(p.knot[l]) $(p.knot[l+1]) $(alph[l+1])\n")
            end
            push!(p.solution.knot, (l, p.knot[l], p.knot[l+1], alph[l+1]))
        end
        if p.legacy_output
            write(fout, "$(p.n*1.0) $(p.knot[p.n]) $(p.knot[p.n+1]) $(alph[p.n] + bet[p.n] * (1-p.knot[p.n]))\n")
            write(fout, "\n")
        end
        push!(p.solution.knot, (p.n, p.knot[p.n], p.knot[p.n+1], alph[p.n] + bet[p.n] * (1 - p.knot[p.n])))

        if p.legacy_output
            write(fout, "  l    alph(l)    bet(l)\n")
        end
        for l = 1:p.n
            if p.legacy_output
                write(fout, "$(l*1.0) $(alph[l]) $(bet[l])\n")
            end
            push!(p.solution.alph_bet, (l, alph[l], bet[l]))
        end
        if p.legacy_output
            write(fout, "\n")
        end

        erre = 0.0
        if p.legacy_output
            write(fout, "  l    x(l)    f(l)\n")
        end
        for l = 1:p.n
            if p.legacy_output
                write(fout, "$(l*1.0) $(x[l]) $(fvec[l])\n")
            end
            push!(p.solution.x_f, (l, x[l], fvec[l]))
            erre += fvec[l]^2
        end
        if p.legacy_output
            write(fout, "\n Norm of derivatives: $(sqrt(erre))\n")
        end
        p.solution.resid = sqrt(erre)

        if p.legacy_output
            close(fout)
        end
    end

    return nothing
end


"""
$(TYPEDSIGNATURES)

Compute CSE for symmetric case defined by `cse_problem`.
"""
function compute_cse(cse_problem::SymmetricAfrprogsCSEProblem, u::Array{Float64})
    @info "Computing: $(cse_problem)"

    # define some arrays
    knot = Vector{Float64}(undef, cse_problem.maxn + 1)
    oldknot = similar(knot)
    yknot = similar(knot, length(knot) + 1)
    x = zeros(Float64, cse_problem.maxn)
    fvec = similar(x)
    alph = similar(x)
    bet = similar(x)

    # parameters initialisation
    # TODO: depend on size on inin
    oldknot[1] = 0.0
    knot[1] = 0.0
    knot[2] = 2.0 / 3.0
    knot[3] = 1.0
    # TODO: move initial guess to cse_problem
    x[1] = log(0.25)
    x[2] = log(0.25)

    # enter a loop that calculates the CSE for different k
    @debug "Entering loop to compute CSE for n=$(cse_problem.inin)..$(cse_problem.maxn)"
    n = cse_problem.inin
    firstiter = true
    solutions = Vector{SymmetricCSESolution}(undef, 0)
    while n <= cse_problem.maxn
        @debug "Loop: n = $n"

        # create Params object for passing extra info to the objective function
        cse_solution = SymmetricCSESolution(problem=cse_problem, n=n)
        prms = SymmetricFunctionParams(n, cse_problem.np, cse_problem.distribution, cse_problem.mc, u, knot, false, cse_solution, cse_problem.legacy_output)

        # solve the system
        # TODO: setup the solver max iters etc the same? fortran values:
        # - max iters = 500
        # - max funevals = 1000
        # - rel error between successive approximations  less than 1e-12
        # - max num iterations 200
        x_n = @view x[begin:n]
        prob = NonlinearProblem(objective_function_symmetric_afrprogs, x_n, prms)
        sol = solve(prob, NewtonRaphson())

        # store some solver info with the solution
        cse_solution.success = SciMLBase.successful_retcode(sol)

        # gather extra details about the solution
        prms.cvrg = true
        fvec_n = @view fvec[begin:n]
        objective_function_symmetric_afrprogs(fvec_n, sol.u, prms)

        # store the solution for this value of n
        push!(solutions, cse_solution)

        # calculate stop criteria C_1 (compare with previous cse)
        if firstiter
            firstiter = false
        else
            csenew = solutions[end].cse."CSE(x)"
            cseold = solutions[end-1].cse."CSE(x)"

            diff = norm(csenew - cseold) / length(csenew)
            solutions[end].c_1 = diff
        end

        # calculate stop criteria C_2 (norm of residual)
        solutions[end].c_2 = norm(sol.resid)

        # log the solution
        @info cse_solution

        # update the parameters before moving to a CSE with a higher k
        yknot[1] = 0.0
        for l = 1:n
            yknot[l+1] = yknot[l] + exp(sol.u[l])
            alph[l] = yknot[l]
            bet[l] = (yknot[l+1] - yknot[l]) / (knot[l+1] - knot[l])
        end

        diff = 0.0
        loc = missing
        for l = 2:n
            oldknot[l] = knot[l]
            aux = abs(bet[l] - bet[l-1])
            if aux > diff
                diff = aux
                loc = l
            end
        end
        oldknot[n+1] = 1.0

        n += 2
        if n <= cse_problem.maxn
            knot[loc+1] = (oldknot[loc-1] + oldknot[loc] + oldknot[loc+1]) / 3.0
            knot[loc] = (oldknot[loc-1] + 2.0 * knot[loc+1]) / 3.0
            knot[n+1] = 1.0

            for l = 1:loc-1
                knot[loc-l+1] = (oldknot[loc-l] + 2.0 * knot[loc+2-l]) / 3.0
            end

            for l = loc+2:n
                knot[l] = (oldknot[l-1] + 2.0 * knot[l-1]) / 3.0
            end

            yknot[1] = 0.0
            for ll in 1:n
                check = true
                l = 0
                ti = knot[ll+1]
                while check
                    l += 1
                    if ti >= oldknot[l] && ti <= oldknot[l+1]
                        yknot[ll+1] = alph[l] + bet[l] * (ti - oldknot[l])
                        check = false
                    end
                end
                x[ll] = log(yknot[ll+1] - yknot[ll])
            end
        end
    end

    return solutions
end
