
"""
$(TYPEDSIGNATURES)

Helper function to call the solver and attach output/stats to the solution
"""
function run_solver(cse_problem::CSEProblem, cse_solution::CSESolution, objective_function::Function, objective_function_params::CSESolverParams, x::AbstractVector{Float64}, fvec::AbstractVector{Float64}, previous_solution::Union{CSESolution,Missing})
    # setting up kwargs for passing to solver
    kwargs = (;)
    if cse_problem.solver_maxiters !== nothing
        kwargs = merge(kwargs, (; maxiters=cse_problem.solver_maxiters))
    end
    if cse_problem.solver_abstol !== nothing
        kwargs = merge(kwargs, (; abstol=cse_problem.solver_abstol))
    end
    if cse_problem.solver_reltol !== nothing
        kwargs = merge(kwargs, (; reltol=cse_problem.solver_reltol))
    end

    # solve the system
    # TODO: setup the solver max iters etc the same? fortran values:
    # - max iters = 500
    # - max funevals = 1000
    # - rel error between successive approximations  less than 1e-12
    # - max num iterations 200
    prob = NonlinearProblem(objective_function, x, objective_function_params)
    sol = solve(prob, cse_problem.solver; kwargs...)

    # store some solver info with the solution
    cse_solution.success = SciMLBase.successful_retcode(sol)
    cse_solution.solver_solution = sol

    if cse_solution.success
        # gather extra details about the solution
        objective_function_params.cvrg = true
        objective_function(fvec, sol.u, objective_function_params)

        # store the solution for this value of n
        # TODO: only push the solutions that succeeded (or make that an option)
        #push!(solutions, cse_solution)

        # TODO: this bit will be different depending on whether a symmetriccsesolution or asym was passed
        # calculate stop criteria C_1 (compare with previous cse)
        if previous_solution !== missing
            csenew = cse_solution.cse."CSE(x)"
            cseold = previous_solution.cse."CSE(x)"

            diff = norm(csenew - cseold) / length(csenew)
            cse_solution.c_1 = diff
        end

        # calculate stop criteria C_2 (norm of residual)
        cse_solution.c_2 = norm(sol.resid)
    else
        @error "Failed to solve CSE => return code is: $(sol.retcode)"
    end

    return sol
end
