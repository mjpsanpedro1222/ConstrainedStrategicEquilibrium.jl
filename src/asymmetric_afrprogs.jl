

# TODO: add stopping criteria tolerances; output files (data, results)
"""
$(TYPEDEF)

Define an asymmetric CSE problem.

Parameters can be passed in as keyword arguments or can be omitted to accept the default values.

$(TYPEDFIELDS)
"""
@kwdef struct AsymmetricAfrprogsCSEProblem <: AsymmetricCSEProblem
    "Random number generator to use during data generation (default rng is seeded with 642867)"
    rng::AbstractRNG = Random.seed!(642867)
    "Number of Monte Carlo steps (default is 10000)"
    mc::Int = 10000
    "Number of players (must be 4 currently)"
    np::Int = 4
    "Distributions to use, which should be a Vector of length `np` (each distribution must be `Beta` currently)"
    distributions::Vector{UnivariateDistribution} = [Beta(3, 3), Beta(3, 3), Beta(5, 3), Beta(5, 3)]
    "Initial value for `n` (default is 16 and must be left as this currently)"
    inin::Int = 16
    "Maximum value for `n` (default is 17)"
    maxn::Int = 17
    "Write txt and csv files with solution info (default is False, most of this info is included in the solution objects that get return from `compute_cse`)"
    legacy_output::Bool = false
    "The solver to use (default is to use `Broyden(; init_jacobian=Val(:true_jacobian), autodiff=AutoFiniteDiff())`)"
    solver::Union{AbstractNonlinearAlgorithm,Nothing} = Broyden(; init_jacobian=Val(:true_jacobian), autodiff=AutoFiniteDiff())
    "Keyword arguments to pass to the solve command, such as `abstol`, `reltol`, `maxiters`, etc. Add `show_trace=Val(true)` to output extra info from the solver."
    solver_kwargs::NamedTuple = (;)
    "Initial guess to pass to the solver, if not provided use a default initial guess (must be length `2 * inin - 1`)"
    solver_initial_guess::Union{Vector{Float64},Nothing} = nothing
end


"""
$(TYPEDSIGNATURES)

Check that `cse_problem` is a valid definition for an asymmetric CSE problem.

This function will throw an error if the problem is not valid, otherwise
will return silently.

# Examples
```jldoctest
julia> prob = AsymmetricAfrprogsCSEProblem();
julia> validate_cse_problem(prob)

julia> prob = AsymmetricAfrprogsCSEProblem(np = 6);
julia> validate_cse_problem(prob)
ERROR: "Only 4 players are supported currently"
[...]

julia> prob = AsymmetricAfrprogsCSEProblem(inin = 4);
julia> validate_cse_problem(prob)
ERROR: "Initial value of n must be 16 currently"
[...]

julia> prob = AsymmetricAfrprogsCSEProblem(maxn = 1);
julia> validate_cse_problem(prob)
ERROR: "Initial value of n cannot be bigger than maximum value of n"
[...]
```
"""
function validate_cse_problem(cse_problem::AsymmetricAfrprogsCSEProblem)
    #    if !(cse_problem.distribution isa Distributions.Beta)
    #        throw("Only Beta distributions are supported currently")
    #    end

    if length(cse_problem.distributions) != cse_problem.np
        throw("Number of distributions must equal the number of players")
    end

    if cse_problem.np != 4
        throw("Only 4 players are supported currently")
    end

    if cse_problem.inin != 16
        throw("Initial value of n must be 16 currently")
    end
    #if cse_problem.inin != cse_problem.maxn
    #    throw("initial value of n must be same as max value currently")
    #end

    if cse_problem.inin > cse_problem.maxn
        throw("Initial value of n cannot be bigger than maximum value of n")
    end

    if cse_problem.solver_initial_guess !== nothing
        expected_length = 2 * cse_problem.inin - 1
        if length(cse_problem.solver_initial_guess) != expected_length
            throw("Solver initial guess must have length $(expected_length) (actual length: $(length(cse_problem.inin)))")
        end
    end
end


"""
$(TYPEDEF)

Structure for passing data to the objective function via the solver interface.
"""
mutable struct AsymmetricFunctionParams <: CSESolverParams
    n::Int64
    np::Int64
    dists::Vector{UnivariateDistribution}
    mc::Int64
    u::Array{Float64,2}
    knot::Array{Float64,2}
    cvrg::Bool
    solution::AsymmetricCSESolution
    legacy_output::Bool
end



"""
$(TYPEDSIGNATURES)

Compute CSE for asymmetric case defined by `cse_problem`.
"""
function compute_cse(cse_problem::AsymmetricAfrprogsCSEProblem, u::Array{Float64})
    @info "Computing: $(cse_problem)"

    # define some arrays
    # TODO: should the sizes be hardcoded or depend on nn? or something else?
    # TODO: can we use static arrays here
    knot = zeros(Float64, 2, cse_problem.maxn + 1)
    oldknot = zeros(Float64, 2, cse_problem.maxn + 1)
    yknot = zeros(Float64, 2, cse_problem.maxn + 2)
    x = zeros(Float64, 2 * cse_problem.maxn - 1)
    fvec = zeros(Float64, 2 * cse_problem.maxn - 1)
    alph = zeros(Float64, 2, cse_problem.maxn + 1)
    bet = zeros(Float64, 2, cse_problem.maxn + 1)

    # initial values
    n = cse_problem.inin
    robo = 15

    # parameters initialisation
    # these were the values provided in the fortran code for n=16
    #    x[1] = -2.84827109173688
    #    x[2] = -2.85209052797506
    #    x[3] = -2.85688597597939
    #    x[4] = -2.86243170560486
    #    x[5] = -2.87005433886236
    #    x[6] = -2.87972857171734
    #    x[7] = -2.89195095174828
    #    x[8] = -2.90986034638539
    #    x[9] = -2.93411212816428
    #    x[10] = -2.96990679761558
    #    x[11] = -3.03212449908756
    #    x[12] = -3.14051739544213
    #    x[13] = -3.33100693581025
    #    x[14] = -3.69184084490141
    #    x[15] = -5.61683383184477
    #    x[16] = -7.41860036191254
    #    x[17] = 4.37486841314759
    #    x[18] = 3.55315379806175
    #    x[19] = 3.09860286531008
    #    x[20] = 2.74440338704673
    #    x[21] = 2.45415125926837
    #    x[22] = 2.19413505997764
    #    x[23] = 1.94448901027584
    #    x[24] = 1.71125414593309
    #    x[25] = 1.47745377147087
    #    x[26] = 1.23478257742541
    #    x[27] = 1.01951187118377
    #    x[28] = 0.846654741348483
    #    x[29] = 0.726435662568677
    #    x[30] = 0.619443991395723
    #    x[31] = 1.12697386952317

    if cse_problem.solver_initial_guess !== nothing && length(cse_problem.solver_initial_guess) == 2 * n - 1
        @debug "Using passed solver initial guess"
        x[begin:2*n-1] .= cse_problem.solver_initial_guess
    elseif cse_problem.solver_initial_guess !== nothing && length(cse_problem.solver_initial_guess) != 2 * n - 1
        @warn "Ignoring passed solver initial guess as it is the wrong length"
    else
        # default values
        @debug "Using default initial guess"
        x[begin:n] .= 2.0 .+ rand(n)
        x[n+1:2*n-1] .= -2.0 .- rand(n - 1)
    end
    @debug "initial guess:" x

    # define the knot array
    tmparr = collect(range(0, stop=1, length=n + 1))
    knot[1, begin:n+1] .= tmparr
    knot[2, begin:n+1] .= tmparr
    @debug "knot 1: $(knot[1, :])"
    @debug "knot 2: $(knot[2, :])"

    # enter a loop that calculates the CSE for different k
    @debug "Entering loop to compute CSE for n=$(cse_problem.inin)..$(cse_problem.maxn)"
    solutions = Vector{AsymmetricCSESolution}(undef, 0)
    previous_solution = missing
    while n <= cse_problem.maxn
        @debug "Loop: n = $n"

        # create Params object for passing extra info to the objective function
        cse_solution = AsymmetricCSESolution(problem=cse_problem, n=n, u=u)
        prms = AsymmetricFunctionParams(n, cse_problem.np, cse_problem.distributions, cse_problem.mc, u, knot, false, cse_solution, cse_problem.legacy_output)

        # solve the system
        x_n = @view x[begin:2*n-1]
        fvec_n = @view fvec[begin:2*n-1]
        #        prob = NonlinearProblem(objective_function_asymmetric_afrprogs, x_n, prms)
        #sol = solve(prob, FastShortcutNonlinearPolyalg(; autodiff=AutoFiniteDiff()))
        #sol = solve(prob, NewtonRaphson(; autodiff=AutoFiniteDiff()))
        #sol = solve(prob, FastShortcutNonlinearPolyalg(; autodiff=AutoFiniteDifferences()))
        #sol = solve(prob)
        #sol = solve(prob, RobustMultiNewton(; autodiff=AutoFiniteDiff()); show_trace=Val(true), trace_level=TraceMinimal(), abstol=1e-4, reltol=1e-4)
        #        sol = solve(prob, RobustMultiNewton(; autodiff=AutoFiniteDiff()); show_trace=Val(true), trace_level=TraceMinimal(), abstol=1e-3, reltol=1e-4)
        sol = run_solver(cse_problem, cse_solution, objective_function_asymmetric_afrprogs, prms, x_n, fvec_n, previous_solution)

        # set previous solution for next step
        previous_solution = cse_solution

        # store the solution for this value of n
        # TODO: only push the solutions that succeeded (or make that an option)
        push!(solutions, cse_solution)

        # log the solution
        @info cse_solution

        # update the parameters before moving to a CSE with a higher k
        yknot[1, 1] = 0.0
        yknot[2, 1] = 0.0
        for l in 1:n-1
            yknot[1, l+1] = yknot[1, l] + exp(sol.u[l])
            alph[1, l] = yknot[1, l]
            bet[1, l] = (yknot[1, l+1] - yknot[1, l]) / (knot[1, l+1] - knot[1, l])
            aux = exp(sol.u[n+l])
            aux = aux / (1 + aux)
            yknot[2, l+1] = yknot[2, l] + (yknot[1, l+1] - yknot[2, l]) * aux
            alph[2, l] = yknot[2, l]
            bet[2, l] = (yknot[2, l+1] - yknot[2, l]) / (knot[2, l+1] - knot[2, l])
        end
        yknot[1, n+1] = yknot[1, n] + exp(sol.u[n])
        alph[1, n] = yknot[1, n]
        bet[1, n] = (yknot[1, n+1] - yknot[1, n]) / (knot[1, n+1] - knot[1, n])
        yknot[2, n+1] = yknot[1, n+1]
        alph[2, n] = yknot[2, n]
        bet[2, n] = (yknot[2, n+1] - yknot[2, n]) / (knot[2, n+1] - knot[2, n])

        n += 1
        if n <= cse_problem.maxn
            oldknot .= knot
            knot[1, n+1-robo] = (knot[1, n-robo] + knot[1, n+1-robo]) / 2.0
            knot[2, n+1-robo] = knot[1, n+1-robo]
            for l in (n+2-robo):(n+1)
                knot[1, l] = oldknot[1, l-1]
                knot[2, l] = knot[1, l]
            end
            robo += 2

            yknot[1, 1] = 0.0
            yknot[2, 1] = 0.0

            for ll in 1:n
                check = true
                l = 0
                ti = knot[1, ll+1]
                while check
                    l += 1
                    if ti >= oldknot[1, l] && ti <= oldknot[1, l+1]
                        yknot[1, ll+1] = alph[1, l] + bet[1, l] * (ti - oldknot[1, l])
                        yknot[2, ll+1] = alph[2, l] + bet[2, l] * (ti - oldknot[1, l])
                        check = false
                    end
                end
                x[ll] = log(yknot[1, ll+1] - yknot[1, ll])
            end

            for l in 1:n-1  # TODO:in their code this was to n, but that is beyond 2n-1???
                aux = (yknot[2, l+1] - yknot[2, l]) / (yknot[1, l+1] - yknot[2, l])
                x[n+l] = log(aux / (1.0 - aux))
            end
        end
    end

    return solutions
end


function objective_function_asymmetric_afrprogs(fvec, x, p::AsymmetricFunctionParams)
    if any(isnan, x)
        @error "NaNs in x at:" findall(isnan, x)
        throw("NaNs in x")
    end

    n = p.n
    n2 = 2 * p.n - 1

    @assert n2 == length(x)
    @assert n2 == length(fvec)

    # note: important to use similar here in case using autodiff they could be of type dual from ForwardDiff
    da = similar(x, length(x) + 2)
    yknot = similar(x, 2, n + 2)
    alph = similar(yknot)
    bet = similar(yknot)

    # distribution
    # NOTE : we assume players 1 and 2 have the same distribution and players 3 and 4 have the same distribution
    #        so we just consider player 1 vs player 3
    dist1 = p.dists[1]
    dist1params = params(dist1)
    dist2 = p.dists[3]
    dist2params = params(dist2)

    # set up the value of the constrained strategy parameters such that the
    # strategy is continuous
    # alph= constant, bet=slope, da(l)= derivative wrt l parameter

    da .= 0.0
    yknot[:, 1] .= 0.0
    for l = 1:n-1
        yknot[1, l+1] = yknot[1, l] + exp(x[l])
        alph[1, l] = yknot[1, l]
        bet[1, l] = (yknot[1, l+1] - yknot[1, l]) / (p.knot[1, l+1] - p.knot[1, l])
        aux = exp(x[n+l])
        aux = aux / (1.0 + aux)
        yknot[2, l+1] = yknot[2, l] + (yknot[1, l+1] - yknot[2, l]) * aux
        alph[2, l] = yknot[2, l]
        bet[2, l] = (yknot[2, l+1] - yknot[2, l]) / (p.knot[2, l+1] - p.knot[2, l])
    end
    yknot[1, n+1] = yknot[1, n] + exp(x[n])
    alph[1, n] = yknot[1, n]
    bet[1, n] = (yknot[1, n+1] - yknot[1, n]) / (p.knot[1, n+1] - p.knot[1, n])

    yknot[2, n+1] = yknot[1, n+1]
    alph[2, n] = yknot[2, n]
    bet[2, n] = (yknot[2, n+1] - yknot[2, n]) / (p.knot[2, n+1] - p.knot[2, n])

    const1 = beta(dist1params...)
    const2 = beta(dist2params...)


    for m = 1:p.mc

        ########################################
        # FOC for Bidder 1
        ########################################

        ti = p.u[m, 1]
        bi = missing
        cumu1 = missing
        dcumu1 = missing
        dbdt1 = missing
        dbdp = missing
        check = true
        l = 0
        while check
            l += 1
            if (ti >= p.knot[1, l]) && (ti <= p.knot[1, l+1])
                bi = alph[1, l] + bet[1, l] * (ti - p.knot[1, l])
                dbdt1 = bet[1, l]
                dbdp = (ti - p.knot[1, l]) / (p.knot[1, l+1] - p.knot[1, l])
                check = false

                cumu1 = cdf(dist1, ti)
                dcumu1 = (ti^(dist1params[1] - 1)) * ((1 - ti)^(dist1params[2] - 1)) / const1
            end
        end

        cumu2 = missing
        dcumu2 = missing
        dbdt2 = missing
        check = true
        ll = 0
        invbi = 1.0
        while check && ll < n
            ll = ll + 1
            if bi >= yknot[2, ll] && bi <= yknot[2, ll+1] && ll <= n
                invbi = (bi - alph[2, ll]) / bet[2, ll] + p.knot[2, ll]
                dbdt2 = bet[2, ll]
                check = false

                cumu2 = cdf(dist2, invbi)
                dcumu2 = (invbi^(dist2params[1] - 1)) * ((1 - invbi)^(dist2params[2] - 1)) / const2
            end
        end

        cumu = cumu1 * cumu2^2
        dcumu = dcumu1 * (cumu2^2) / dbdt1 + 2 * dcumu2 * cumu1 * cumu2 / dbdt2
        da[l] = da[l] + dbdp * ((ti - bi) * dcumu - cumu)

        ########################################
        # FOC for Bidder 4
        ########################################

        ti = p.u[m, 4]
        bi = missing
        check = true
        l = 0
        while check
            l += 1
            if ti >= p.knot[2, l] && ti <= p.knot[2, l+1]
                bi = alph[2, l] + bet[2, l] * (ti - p.knot[2, l])
                dbdt1 = bet[2, l]
                dbdp = (ti - p.knot[2, l]) / (p.knot[2, l+1] - p.knot[2, l])
                dbdp = dbdp * (p.knot[1, l+1] - p.knot[2, l])
                check = false

                cumu1 = cdf(dist2, ti)
                dcumu1 = (ti^(dist2params[1] - 1)) * ((1 - ti)^(dist2params[2] - 1)) / const2
            end
        end

        check = true
        ll = 0
        invbi = 1.0
        while check && ll < n
            ll += 1
            if bi >= yknot[1, ll] && bi <= yknot[1, ll+1] && ll <= n
                invbi = (bi - alph[1, ll]) / bet[1, ll] + p.knot[1, ll]
                dbdt2 = bet[1, ll]
                check = false

                cumu2 = cdf(dist1, invbi)
                dcumu2 = (invbi^(dist1params[1] - 1)) * ((1 - invbi)^(dist1params[2] - 1)) / const1
            end
        end

        cumu = cumu1 * cumu2^2
        dcumu = dcumu1 * (cumu2^2) / dbdt1 + 2 * dcumu2 * cumu1 * cumu2 / dbdt2
        da[n+l] += dbdp * ((ti - bi) * dcumu - cumu)
    end

    for l = 1:p.n-1
        fvec[l] = da[l] * yknot[1, l+1] / p.mc
        aux = exp(x[n+l])
        aux = aux / (1.0 + aux)^2
        fvec[n+l] = da[n+l] * aux / p.mc
    end
    fvec[p.n] = da[p.n] * yknot[1, n+1] / p.mc
    if any(isnan, fvec)
        throw("nans in fvec")
    end

    if p.cvrg
        if p.legacy_output
            fout = open("asym-result-n-$(p.n).txt", "w")
            fcsv = open("asym-bids-private-values-n-$(p.n).csv", "w")

            for m = 1:2*n-1
                write(fout, "x($(m)) = $(x[m])\n")
            end

            for m = 1:n+1
                write(fout, "knot(1, $(m)) = $(p.knot[1, m])\n")
            end
            write(fout, "\n")

            write(fout, "t    CSE(t)\n")
            write(fcsv, "t,CSE(t) 1,CSE(t) 2\n")
        end

        for m = 1:101
            ti = (m - 1.0) / 100.0
            bi = missing
            check = true
            l = 0
            while check
                l += 1
                if (ti >= p.knot[1, l] && ti <= p.knot[1, l+1])
                    bi = alph[1, l] + bet[1, l] * (ti - p.knot[1, l])
                    check = false
                end
            end

            check = true
            bbi = missing
            l = 0
            while check
                l += 1
                if (ti >= p.knot[2, l] && ti <= p.knot[2, l+1])
                    bbi = alph[2, l] + bet[2, l] * (ti - p.knot[2, l])
                    check = false
                end
            end

            if p.legacy_output
                write(fout, "$(ti)    $(bi)   $(bbi)\n")
                write(fcsv, "$ti,$bi,$bbi\n")
            end
            push!(p.solution.cse, (ti, bi, bbi))
        end

        if p.legacy_output
            write(fout, "\n")
            close(fcsv)
        end

        # TODO: store these on the solution too
        if p.legacy_output
            write(fout, "PLAYERS 1 and 2\n")
            write(fout, "  l    knot(1,l-1)    knot(1,l)    CSE[knot(1,l)]\n")
            for l = 1:p.n-1
                write(fout, "$(l*1.0) $(p.knot[1, l]) $(p.knot[1, l+1]) $(alph[1, l+1])\n")
            end
            write(fout, "$(p.n*1.0) $(p.knot[1, p.n]) $(p.knot[1, p.n+1]) $(alph[1, p.n] + bet[1, p.n] * (1-p.knot[1, p.n]))\n")
            write(fout, "\n")

            write(fout, "  l    alph(1,l)    bet(1,l)\n")
            for l = 1:p.n
                write(fout, "$(l*1.0) $(alph[1, l]) $(bet[1, l])\n")
            end
            write(fout, "\n")

            erre = 0.0
            write(fout, "  l    x(l)    f(l)\n")
            for l = 1:p.n
                write(fout, "$(l*1.0) $(x[l]) $(fvec[l])\n")
                erre += fvec[l]^2
            end
            write(fout, "Norm of derivatives: $(sqrt(erre))\n")

            write(fout, "\n\n")
            write(fout, "PLAYERS 3 and 4\n")
            write(fout, "  l    knot(2,l-1)    knot(2,l)    CSE[knot(2,l)]\n")
            for l = 1:p.n-1
                write(fout, "$(l*1.0) $(p.knot[2, l]) $(p.knot[2, l+1]) $(alph[2, l+1])\n")
            end
            write(fout, "$(p.n*1.0) $(p.knot[2, p.n]) $(p.knot[2, p.n+1]) $(alph[2, p.n] + bet[2, p.n] * (1-p.knot[2, p.n]))\n")
            write(fout, "\n")

            write(fout, "  l    alph(2,l)    bet(2,l)\n")
            for l = 1:p.n
                write(fout, "$(l*1.0) $(alph[2, l]) $(bet[2, l])\n")
            end
            write(fout, "\n")

            erre = 0.0
            write(fout, "  l    x(l)    f(l)\n")
            for l = 1:p.n-1
                write(fout, "$(l*1.0) $(x[p.n+l]) $(fvec[p.n+l])\n")
                erre += fvec[p.n+l]^2
            end
            write(fout, "Norm of derivatives: $(sqrt(erre))\n")

            close(fout)
        end
    end

    return nothing
end
