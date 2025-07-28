
"""
$(TYPEDSIGNATURES)

Compute the given CSE problem.

This function will validate the problem, generate the data and the compute the CSE, returning the vector of solutions.
"""
function compute_cse(cse_problem::CSEProblem)
    # validate the problem
    validate_cse_problem(cse_problem)

    # generate the data
    u = generate_cse_data(cse_problem)

    # call the underlying implementation
    result = compute_cse(cse_problem, u)

    return result
end
