

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
    "Number of players (must be 2 or 4 currently)"
    np::Int = 2
    "Distributions to use, which should be a Vector of length `np` (valid distributions are Beta and Uniform currently)"
    distributions::Vector{UnivariateDistribution} = [Beta(3, 3), Beta(5, 3)]
    "Initial value for `n` (default is 16)"
    inin::Int = 16
    "Maximum value for `n` (default is 17)"
    maxn::Int = 17
    "Knot refinement strategy. Can be `:steepest_slope`, `:highest_curvature`, `:even_spacing`, or `:double_knot`. (default is `:steepest_slope`)"
    knot_refinement_strategy::Symbol = :steepest_slope
    "Write txt and csv files with solution info (default is False, most of this info is included in the solution objects that get return from `compute_cse`)"
    legacy_output::Bool = false
    "The solver to use (default is to use `Broyden(; init_jacobian=Val(:true_jacobian), autodiff=AutoFiniteDiff())`)"
    solver::Union{AbstractNonlinearAlgorithm,Nothing} = Broyden(; init_jacobian=Val(:true_jacobian), autodiff=AutoFiniteDiff())
    "Keyword arguments to pass to the solve command, such as `abstol`, `reltol`, `maxiters`, etc. Add `show_trace=Val(true)` to output extra info from the solver."
    solver_kwargs::NamedTuple = (;)
    "Initial guess to pass to the solver, if not provided use a default initial guess (must be length `2 * inin - 1`)"
    solver_initial_guess::Union{Vector{Float64},Nothing} = nothing
    "Initial knot positions to use, if not provided use a default initial guess (must be length `inin + 1`)"
    initial_knots::Union{Vector{Float64},Nothing} = nothing
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
ERROR: "Only 2 or 4 players are supported currently"
[...]

julia> prob = AsymmetricAfrprogsCSEProblem(maxn = 1);
julia> validate_cse_problem(prob)
ERROR: "Initial value of n cannot be bigger than maximum value of n"
[...]
```
"""
function validate_cse_problem(cse_problem::AsymmetricAfrprogsCSEProblem)
    if !all(x -> isa(x, Distributions.Beta) || isa(x, Distributions.Uniform), cse_problem.distributions)
        throw("Only Beta or Uniform distributions are supported currently")
    end

    if length(cse_problem.distributions) != cse_problem.np
        if (length(cse_problem.distributions) == 2) && (cse_problem.np == 4)
            @info "Using $(cse_problem.distributions[1]) for bidders 1 and 2; $(cse_problem.distributions[2]) for bidders 3 and 4."
            resize!(cse_problem.distributions, 4)
            cse_problem.distributions[3] = cse_problem.distributions[2]
            cse_problem.distributions[4] = cse_problem.distributions[2]
            cse_problem.distributions[2] = cse_problem.distributions[1]
        else
            throw("Number of distributions must equal the number of players")
        end
    end

    if length(cse_problem.distributions) == 4
        if (cse_problem.distributions[1] != cse_problem.distributions[2]) || (cse_problem.distributions[3] != cse_problem.distributions[4])
            throw("Bidders 1 and 2 must have the same distribution and the same goes for bidders 3 and 4")
        end
    end

    if (cse_problem.np != 4) && (cse_problem.np != 2)
        throw("Only 2 or 4 players are supported currently")
    end

    if cse_problem.inin < 1
        throw("Initial value of n must greater than 0")
    end

    if cse_problem.inin > cse_problem.maxn
        throw("Initial value of n cannot be bigger than maximum value of n")
    end

    if cse_problem.knot_refinement_strategy ∉ [:steepest_slope, :highest_curvature, :even_spacing, :double_knot]
        throw("`knot_refinement_strategy` must be one of `:steepest_slope`, `:highest_curvature`, `:even_spacing`, or `:double_knot`")
    end

    if cse_problem.solver_initial_guess !== nothing
        expected_length = 2 * cse_problem.inin - 1
        if length(cse_problem.solver_initial_guess) != expected_length
            throw("Solver initial guess must have length $(expected_length) (actual length: $(length(cse_problem.solver_initial_guess)))")
        end
    end

    if cse_problem.initial_knots !== nothing
        expected_length = cse_problem.inin + 1
        if length(cse_problem.initial_knots) != expected_length
            throw("Initial knots must have length $(expected_length) (actual length: $(length(cse_problem.initial_knots)))")
        end
        if !issorted(cse_problem.initial_knots) || cse_problem.initial_knots[begin] != 0.0 || cse_problem.initial_knots[end] != 1.0
            throw("Initial knots must be sorted and range from 0.0 to 1.0")
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
    knot = zeros(Float64, 2, cse_problem.maxn + 1)
    oldknot = zeros(Float64, 2, cse_problem.maxn + 1)
    yknot = zeros(Float64, 2, cse_problem.maxn + 2)
    x = zeros(Float64, 2 * cse_problem.maxn - 1)
    fvec = zeros(Float64, 2 * cse_problem.maxn - 1)
    alph = zeros(Float64, 2, cse_problem.maxn + 1)
    bet = zeros(Float64, 2, cse_problem.maxn + 1)

    # initial values
    n = cse_problem.inin

    # parameters initialisation
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

    # define the knot array
    if cse_problem.initial_knots !== nothing
        knot[1, begin:n+1] .= cse_problem.initial_knots
    else
        knot[1, begin:n+1] .= range(0, stop=1, length=n + 1)
    end
    knot[2, begin:n+1] .= knot[1, begin:n+1]

    # enter a loop that calculates the CSE for different k
    @debug "Entering loop to compute CSE for n=$(cse_problem.inin)..$(cse_problem.maxn)"
    solutions = Vector{AsymmetricCSESolution}(undef, 0)
    previous_solution = missing
    while n <= cse_problem.maxn
        @debug "Loop: n = $n"
        @debug "Initial guess:" x
        @debug "Knot 1:" knot[1, :]
        @debug "Knot 2:" knot[2, :]

        # create Params object for passing extra info to the objective function
        cse_solution = AsymmetricCSESolution(problem=cse_problem, n=n, u=u)
        prms = AsymmetricFunctionParams(n, cse_problem.np, cse_problem.distributions, cse_problem.mc, u, knot, false, cse_solution, cse_problem.legacy_output)

        # solve the system
        x_n = @view x[begin:2*n-1]
        fvec_n = @view fvec[begin:2*n-1]
        sol = run_solver(cse_problem, cse_solution, objective_function_asymmetric_afrprogs, prms, x_n, fvec_n, previous_solution)

        # set previous solution for next step
        previous_solution = cse_solution

        # store the solution for this value of n
        push!(solutions, cse_solution)

        # log the solution
        @info cse_solution

        # break the loop if failed
        if !cse_solution.success
            @error "Exiting compute_cse due to solve failed"
            break
        end

        newn = (cse_problem.knot_refinement_strategy == :double_knot) ? n + 2 : n + 1
        if newn <= cse_problem.maxn
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

            if cse_problem.knot_refinement_strategy == :double_knot
                # Use symmetric approach: add 2 knots and use 1/3, 2/3 positioning
                # find location with highest curvature (like symmetric)
                diff = 0.0
                loc = 2
                for l = 2:n
                    aux = abs(bet[1, l] - bet[1, l-1]) + abs(bet[2, l] - bet[2, l-1])
                    if aux > diff
                        diff = aux
                        loc = l
                    end
                end

                n += 2
                oldknot .= knot

                # Apply symmetric-style knot positioning
                knot[1, loc+1] = (oldknot[1, loc-1] + oldknot[1, loc] + oldknot[1, loc+1]) / 3.0
                knot[1, loc] = (oldknot[1, loc-1] + 2.0 * knot[1, loc+1]) / 3.0
                knot[1, n+1] = 1.0

                # Update knots before the split location
                for l = 1:loc-1
                    knot[1, loc-l+1] = (oldknot[1, loc-l] + 2.0 * knot[1, loc+2-l]) / 3.0
                end

                # Update knots after the split location
                for l = loc+2:n
                    knot[1, l] = (oldknot[1, l-1] + 2.0 * knot[1, l-1]) / 3.0
                end

                # Copy to second player (maintaining assumption they're the same)
                knot[2, 1:n+1] .= knot[1, 1:n+1]
            else
                # find interval to split for next iteration
                split_idx = -1
                if cse_problem.knot_refinement_strategy == :steepest_slope
                    # split interval with largest derivative (steepest slope)
                    slopes = abs.(bet[1, 1:n]) + abs.(bet[2, 1:n])
                    ~, split_idx = findmax(slopes)
                elseif cse_problem.knot_refinement_strategy == :highest_curvature
                    # split interval with largest change in derivative (highest curvature)
                    if n > 1
                        bet_diff1 = abs.(bet[1, 2:n] - bet[1, 1:n-1])
                        bet_diff2 = abs.(bet[2, 2:n] - bet[2, 1:n-1])
                        total_diffs = bet_diff1 + bet_diff2
                        ~, idx = findmax(total_diffs)

                        # between the two intervals that make the corner, pick the one with the steeper slope
                        slopes = abs.(bet[1, 1:n]) + abs.(bet[2, 1:n])
                        if slopes[idx] >= slopes[idx+1]
                            split_idx = idx
                        else
                            split_idx = idx + 1
                        end
                    else
                        # cannot compute curvature with one interval, fallback to steepest slope
                        slopes = abs.(bet[1, 1:n]) + abs.(bet[2, 1:n])
                        ~, split_idx = findmax(slopes)
                    end
                end

                n += 1
                oldknot .= knot

                if cse_problem.knot_refinement_strategy == :even_spacing
                    @debug "Regenerating evenly spaced knots for n=$(n)"
                    knot[1, begin:n+1] .= range(0, stop=1, length=n + 1)
                    knot[2, begin:n+1] .= knot[1, begin:n+1]
                else
                    # insert new knot point by splitting the interval with the largest change
                    @debug "Inserting new knot point after: $(split_idx)"
                    new_knot_val = (oldknot[1, split_idx] + oldknot[1, split_idx+1]) / 2.0
                    knot[1, split_idx+1] = new_knot_val
                    knot[1, split_idx+2:n+1] .= oldknot[1, split_idx+1:n]
                    # NOTE: there are assumptions elsewhere in the code about the knot points being the same for both players
                    knot[2, 1:n+1] .= knot[1, 1:n+1]
                end
            end

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
                if aux == 1.0
                    aux = 1 - 1e-15
                end
                x[n+l] = log(aux / (1.0 - aux))
                #@debug "calc next x" n+l x[n+l] aux yknot[2, l+1] yknot[2, l] yknot[1, l+1]
            end
        else
            n += 1
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
    da .= 0
    yknot = similar(x, 2, n + 2)
    yknot .= 0
    alph = similar(yknot)
    alph .= 0
    bet = similar(yknot)
    bet .= 0

    # distribution
    # NOTE : in the 4 player case, we assume players 1 and 2 have the same distribution and
    #        players 3 and 4 have the same distribution so we just consider player 1 vs player 3
    dist1 = p.dists[1]
    if p.np == 2
        dist2 = p.dists[2]
    else
        dist2 = p.dists[3]
    end

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

    if any(isnan, yknot)
        @error "NaNs in yknot at:" findall(isnan, yknot) x
        throw("NaNs in yknot")
    end

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
                dcumu1 = pdf(dist1, ti)
            end
        end
        if cumu1 === missing
            @error "Failed to find interval for ti in FOC 1" ti p.knot
            throw("Failed to find interval for ti in FOC 1")
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
                dcumu2 = pdf(dist2, invbi)
            end
        end
        if cumu2 === missing
            @error "Failed to find interval for bi in FOC 1" bi yknot
            throw("Failed to find interval for bi in FOC 1")
        end

        # Number of other bidders of same type (type 1)
        n_same_type = (p.np == 2) ? 0 : 1  # 0 for 2-player, 1 for 4-player
        # Number of bidders of other type (type 2) 
        n_other_type = (p.np == 2) ? 1 : 2  # 1 for 2-player, 2 for 4-player

        if n_same_type == 0
            cumu = cumu2^n_other_type
            dcumu = dcumu1 * cumu2^n_other_type / dbdt1 + n_other_type * dcumu2 * cumu2^(n_other_type-1) / dbdt2
        else
            cumu = cumu1^n_same_type * cumu2^n_other_type
            dcumu = n_same_type * dcumu1 * cumu1^(n_same_type-1) * cumu2^n_other_type / dbdt1 + 
                    n_other_type * dcumu2 * cumu1^n_same_type * cumu2^(n_other_type-1) / dbdt2
        end
        da[l] = da[l] + dbdp * ((ti - bi) * dcumu - cumu)

        ########################################
        # FOC for Bidder 4
        ########################################

        if p.np == 2
            ti = p.u[m, 2]
        else
            ti = p.u[m, 4]
        end
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
                dcumu1 = pdf(dist2, ti)
            end
        end
        if bi === missing
            @error "Failed to find interval for ti in FOC 2" ti p.knot
            throw("Failed to find interval for ti in FOC 2")
        end

        check = true
        ll = 0
        invbi = 1.0
        cumu2 = missing
        while check && ll < n
            ll += 1
            if bi >= yknot[1, ll] && bi <= yknot[1, ll+1] && ll <= n
                invbi = (bi - alph[1, ll]) / bet[1, ll] + p.knot[1, ll]
                dbdt2 = bet[1, ll]
                check = false

                cumu2 = cdf(dist1, invbi)
                dcumu2 = pdf(dist1, invbi)
            end
        end
        if cumu2 === missing
            @error "Failed to find interval for bi in FOC 2" bi yknot
            throw("Failed to find interval for bi in FOC 2")
        end

        # Number of other bidders of same type (type 2)
        n_same_type = (p.np == 2) ? 0 : 1  # 0 for 2-player, 1 for 4-player
        # Number of bidders of other type (type 1)
        n_other_type = (p.np == 2) ? 1 : 2  # 1 for 2-player, 2 for 4-player

        if n_same_type == 0
            cumu = cumu2^n_other_type
            dcumu = dcumu1 * cumu2^n_other_type / dbdt1 + n_other_type * dcumu2 * cumu2^(n_other_type-1) / dbdt2
        else
            cumu = cumu1^n_same_type * cumu2^n_other_type
            dcumu = n_same_type * dcumu1 * cumu1^(n_same_type-1) * cumu2^n_other_type / dbdt1 + 
                    n_other_type * dcumu2 * cumu1^n_same_type * cumu2^(n_other_type-1) / dbdt2
        end
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

            write(fout, "PLAYERS 1 and 2\n")
            write(fout, "  l    knot(1,l-1)    knot(1,l)    CSE[knot(1,l)]\n")
        end
        for l = 1:p.n-1
            push!(p.solution.knot[:bidder1], (l, p.knot[1, l], p.knot[1, l+1], alph[1, l+1]))
            if p.legacy_output
                write(fout, "$(l*1.0) $(p.knot[1, l]) $(p.knot[1, l+1]) $(alph[1, l+1])\n")
            end
        end
        push!(p.solution.knot[:bidder1], (p.n, p.knot[1, p.n], p.knot[1, p.n+1], alph[1, p.n] + bet[1, p.n] * (1 - p.knot[1, p.n])))
        if p.legacy_output
            write(fout, "$(p.n*1.0) $(p.knot[1, p.n]) $(p.knot[1, p.n+1]) $(alph[1, p.n] + bet[1, p.n] * (1-p.knot[1, p.n]))\n")
            write(fout, "\n")

            write(fout, "  l    alph(1,l)    bet(1,l)\n")
        end
        for l = 1:p.n
            push!(p.solution.alph_bet[:bidder1], (l, alph[1, l], bet[1, l]))
            if p.legacy_output
                write(fout, "$(l*1.0) $(alph[1, l]) $(bet[1, l])\n")
            end
        end
        if p.legacy_output
            write(fout, "\n")

            write(fout, "  l    x(l)    f(l)\n")
        end
        erre = 0.0
        for l = 1:p.n
            push!(p.solution.x_f[:bidder1], (l, x[l], fvec[l]))
            if p.legacy_output
                write(fout, "$(l*1.0) $(x[l]) $(fvec[l])\n")
            end
            erre += fvec[l]^2
        end
        p.solution.norm_derivatives[:bidder1] = sqrt(erre)
        if p.legacy_output
            write(fout, "Norm of derivatives: $(sqrt(erre))\n")

            write(fout, "\n\n")
            write(fout, "PLAYERS 3 and 4\n")
            write(fout, "  l    knot(2,l-1)    knot(2,l)    CSE[knot(2,l)]\n")
        end
        for l = 1:p.n-1
            push!(p.solution.knot[:bidder2], (l, p.knot[2, l], p.knot[2, l+1], alph[2, l+1]))
            if p.legacy_output
                write(fout, "$(l*1.0) $(p.knot[2, l]) $(p.knot[2, l+1]) $(alph[2, l+1])\n")
            end
        end
        push!(p.solution.knot[:bidder2], (p.n, p.knot[2, p.n], p.knot[2, p.n+1], alph[2, p.n] + bet[2, p.n] * (1 - p.knot[2, p.n])))
        if p.legacy_output
            write(fout, "$(p.n*1.0) $(p.knot[2, p.n]) $(p.knot[2, p.n+1]) $(alph[2, p.n] + bet[2, p.n] * (1-p.knot[2, p.n]))\n")
            write(fout, "\n")

            write(fout, "  l    alph(2,l)    bet(2,l)\n")
        end
        for l = 1:p.n
            push!(p.solution.alph_bet[:bidder2], (l, alph[2, l], bet[2, l]))
            if p.legacy_output
                write(fout, "$(l*1.0) $(alph[2, l]) $(bet[2, l])\n")
            end
        end
        if p.legacy_output
            write(fout, "\n")

            write(fout, "  l    x(l)    f(l)\n")
        end
        erre = 0.0
        for l = 1:p.n-1
            push!(p.solution.x_f[:bidder2], (l, x[p.n+l], fvec[p.n+l]))
            if p.legacy_output
                write(fout, "$(l*1.0) $(x[p.n+l]) $(fvec[p.n+l])\n")
            end
            erre += fvec[p.n+l]^2
        end
        p.solution.norm_derivatives[:bidder2] = sqrt(erre)
        if p.legacy_output
            write(fout, "Norm of derivatives: $(sqrt(erre))\n")

            close(fout)
        end
    end

    return nothing
end
