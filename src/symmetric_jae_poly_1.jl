
"""
$(TYPEDEF)

The jae_poly_1 symmetric CSE problem from Computer_Code_CSE.

Parameters can be passed in as keyword arguments or can be omitted to accept the default values.

$(TYPEDFIELDS)

# Examples
```jldoctest
julia> prob = SymmetricJaePoly1CSEProblem()
SymmetricJaePoly1CSEProblem(np=4, mc=10000, n=2..16, Distributions.Beta{Float64}(α=3.0, β=3.0))

julia> prob = SymmetricJaePoly1CSEProblem(mc = 1000, maxn = 12, distribution = Beta(3, 4))
SymmetricJaePoly1CSEProblem(np=4, mc=1000, n=2..12, Distributions.Beta{Float64}(α=3.0, β=4.0))
```
"""
@kwdef struct SymmetricJaePoly1CSEProblem <: SymmetricCSEProblem
    "Random number generator to use during data generation (default rng is seeded with 642867)"
    rng::AbstractRNG = Random.seed!(642867)
    "Number of Monte Carlo steps (default is 10000)"
    mc::Int = 10000
    "Number of players (default is 4)"
    np::Int = 4
    "Distribution to use (must be `Kumaraswamy` currently; default is `Kumaraswamy(2.5, 3.5)`)"
    distribution::UnivariateDistribution = Kumaraswamy(2.5, 3.5)
    "Initial value for n (default is 1 and must be left as this currently)"
    inin::Int = 1
    "Maximum value for n (default is 5)"
    maxn::Int = 5
    "Write txt and csv files with solution info (default is false)"
    legacy_output::Bool = false
    "The solver to use (default is to use the default set by NonlinearSolve.jl)"
    solver::Union{AbstractNonlinearAlgorithm,Nothing} = nothing
    "Keyword arguments to pass to the solve command, such as abstol, reltol, maxiters, etc."
    solver_kwargs::NamedTuple = (;)
end


"""
$(TYPEDSIGNATURES)

Check that `cse_problem` is a valid definition for a CSE problem.

This function will throw an error if the problem is not valid, otherwise
will return silently.

# Examples
```jldoctest
julia> prob = SymmetricJaePoly1CSEProblem();
julia> validate_cse_problem(prob)

julia> prob = SymmetricJaePoly1CSEProblem(distribution = Normal());
julia> validate_cse_problem(prob)
ERROR: "Only Kumaraswamy distributions are supported currently"
[...]

julia> prob = SymmetricJaePoly1CSEProblem(inin = 4);
julia> validate_cse_problem(prob)
ERROR: "Initial value of n must be 1 currently"
[...]

julia> prob = SymmetricJaePoly1CSEProblem(maxn = 0);
julia> validate_cse_problem(prob)
ERROR: "Initial value of n cannot be bigger than maximum value of n"
[...]
```
"""
function validate_cse_problem(cse_problem::SymmetricJaePoly1CSEProblem)
    if !(cse_problem.distribution isa Distributions.Kumaraswamy)
        throw("Only Kumaraswamy distributions are supported currently")
    end

    if cse_problem.np < 2
        throw("Not enough players")
    end

    if cse_problem.inin != 1
        throw("Initial value of n must be 1 currently")
    end

    if cse_problem.inin > cse_problem.maxn
        throw("Initial value of n cannot be bigger than maximum value of n")
    end
end


mutable struct PolyParams <: CSESolverParams
    n::Int64
    np::Int64
    dist::UnivariateDistribution
    mc::Int64
    u::Array{Float64,2}
    knot::Float64
    cvrg::Bool
    solution::SymmetricCSESolution
    legacy_output::Bool
end


"""
$(TYPEDSIGNATURES)

Objective function for the symmetric jae_poly_1 case.
"""
function objective_function_symmetric_jaepoly1(fvec, x, p::PolyParams)
    # note: important to use similar here in case using autodiff they could be of type dual from ForwardDiff
    da = similar(x)

    da .= 0.0
    a = 0.0
    b = 0.0
    for l in 1:p.n
        b += l * x[l] * p.knot^(l - 1)
        a += x[l] * p.knot^l
    end

    # parameters of the distribution
    distprms = params(p.dist)

    for m = 1:p.mc
        ti = p.u[m, 1]  # just for the first player
        if ti <= p.knot
            cumu1 = 1.0 - (1.0 - ti^distprms[1])^distprms[2]
            dcumu1 = distprms[1] * distprms[2] * (ti^(distprms[1] - 1.0)) * (1.0 - ti^distprms[1])^(distprms[2] - 1.0)
            cumu = cumu1^(p.np - 1)
            dcumu = (p.np - 1.0) * dcumu1 * cumu1^(p.np - 2.0)

            bi = 0.0
            dbdti = 0.0
            for l in 1:p.n
                dbdti += x[l] * l * ti^(l - 1)
                bi += x[l] * ti^l
            end

            if dbdti < 0
                dbdti = abs(dbdti)
            end

            for l in 1:p.n
                da[l] += (ti^l) * ((ti - bi) * dcumu / dbdti - cumu)
            end
        end
    end

    fvec .= da ./ p.mc

    # if cvrg flag is set, output results
    if p.cvrg
        if p.legacy_output
            fout = open("sym-result-n-$(p.n).txt", "w")
            fcsv = open("sym-bids-private-values-n-$(p.n).csv", "w")
            write(fout, "  t    CSE(t)\n")
            write(fcsv, "t,CSE(t)\n")
        end

        mse_sum = 0.0
        count_valid = 0
        for m = 1:101
            ti = (m - 1.0) / 100.0
            bi = 0.0

            if ti <= p.knot
                for l = 1:p.n
                    bi += x[l] * ti^l
                end
            else
                bi = a + b * (ti - p.knot)
            end

            true_bne = compute_bne(ti, p.dist, p.np)

            # If the first element's result is NaN, replace it with 0.0
            if m == 1 && isnan(true_bne)
                true_bne = 0.0
            end

            # Compute squared error if bi is not missing
            if !ismissing(bi)
                mse_sum += (bi - true_bne)^2
                count_valid += 1
            end

            if p.legacy_output
                write(fout, "$ti    $bi\n")
                write(fcsv, "$ti,$bi\n")
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

Compute CSE for "jae_poly_1" symmetric case defined by `cse_problem`.
"""
function compute_cse(cse_problem::SymmetricJaePoly1CSEProblem, u::Array{Float64})
    @info "Computing: $(cse_problem)"

    # define some arrays
    x = zeros(Float64, cse_problem.maxn)
    fvec = zeros(Float64, cse_problem.maxn)

    # parameters initialisation
    # TODO: move initial guess to cse_problem
    x[1] = 0.5
    knot = 0.95

    # enter a loop that calculates the CSE for different k
    @debug "Entering loop to compute CSE for n=$(cse_problem.inin)..$(cse_problem.maxn)"
    n = cse_problem.inin
    solutions = Vector{SymmetricCSESolution}(undef, 0)
    previous_solution = missing
    while n <= cse_problem.maxn
        @debug "Loop: n = $n"

        # create Params object for passing extra info to the objective function
        cse_solution = SymmetricCSESolution(problem=cse_problem, n=n, u=u)
        prms = PolyParams(n, cse_problem.np, cse_problem.distribution, cse_problem.mc, u, knot, false, cse_solution, cse_problem.legacy_output)

        # solve the system
        x_n = @view x[begin:n]
        fvec_n = @view fvec[begin:n]
        sol = run_solver(cse_problem, cse_solution, objective_function_symmetric_jaepoly1, prms, x_n, fvec_n, previous_solution)

        # set previous solution to latest ready for next step
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

        n += 1
        if n <= cse_problem.maxn
            x_n .= sol.u
            x[n] = 0.0
        end
    end

    return solutions
end
