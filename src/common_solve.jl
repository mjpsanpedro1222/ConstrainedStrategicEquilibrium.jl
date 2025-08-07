
"""
$(TYPEDSIGNATURES)

Helper function to call the solver and attach output/stats to the solution
"""
function run_solver(cse_problem::CSEProblem, cse_solution::CSESolution, objective_function::Function, objective_function_params::CSESolverParams, x::AbstractVector{Float64}, fvec::AbstractVector{Float64}, previous_solution::Union{CSESolution,Missing})
    # solve the system
    prob = NonlinearProblem(objective_function, x, objective_function_params)

    sol = nothing
    try
        sol = solve(prob, cse_problem.solver; cse_problem.solver_kwargs...)
    catch e
        @error "Caught error while running the solver" e
        cse_solution.success = false
        cse_solution.solver_solution = nothing
    else
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
                elseif isa(cse_solution, AsymmetricCSESolution)
                    # player 1
                    csenew = cse_solution.cse."CSE(x) 1"
                    cseold = previous_solution.cse."CSE(x) 1"
                    diff = norm(csenew - cseold) / length(csenew)
                    cse_solution.c_1[:bidder1] = diff

                    # player 2
                    csenew = cse_solution.cse."CSE(x) 2"
                    cseold = previous_solution.cse."CSE(x) 2"
                    diff = norm(csenew - cseold) / length(csenew)
                    cse_solution.c_1[:bidder2] = diff
                else
                    @warn "C_1 not implemented for $(typeof(cse_solution))"
                end
            end

            # calculate stop criteria C_2 (norm of residual)
            cse_solution.c_2 = norm(sol.resid)
        else
            @error "Failed to solve CSE => return code is: $(sol.retcode)"
        end
    end

    return sol
end
