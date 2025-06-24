
"""
$(TYPEDSIGNATURES)

Helper function to call the solver and attach output/stats to the solution
"""
function run_solver(cse_problem::CSEProblem, cse_solution::CSESolution, objective_function::Function, objective_function_params::CSESolverParams, x::AbstractVector{Float64}, fvec::AbstractVector{Float64}, previous_solution::Union{CSESolution,Missing})
    # solve the system
    prob = NonlinearProblem(objective_function, x, objective_function_params)
    sol = solve(prob, cse_problem.solver; cse_problem.solver_kwargs...)

    # store some solver info with the solution
    cse_solution.success = SciMLBase.successful_retcode(sol)
    cse_solution.solver_solution = sol

    if cse_solution.success
        # gather extra details about the solution
        objective_function_params.cvrg = true
        objective_function(fvec, sol.u, objective_function_params)

        # calculate stop criteria C_1 (compare with previous cse)
        if previous_solution !== missing
            if isa(cse_solution, SymmetricCSESolution)
                csenew = cse_solution.cse."CSE(x)"
                cseold = previous_solution.cse."CSE(x)"

                diff = norm(csenew - cseold) / length(csenew)
                cse_solution.c_1 = diff
            elseif isa(cse_solution, AsymmetricCSEProblem)
                # TODO: just for group 1 currently, add for all players eventually?
                csenew = cse_solution.cse."CSE(x) group 1"
                cseold = previous_solution.cse."CSE(x) group 1"

                diff = norm(csenew - cseold) / length(csenew)
                cse_solution.c_1 = diff
            else
                @warn "C_1 not implemented for $(typeof(cse_solution))"
            end
        end

        # calculate stop criteria C_2 (norm of residual)
        cse_solution.c_2 = norm(sol.resid)
    else
        @error "Failed to solve CSE => return code is: $(sol.retcode)"
    end

    return sol
end
