
function compute_cse(cse_problem::SymmetricCSEProblem)
    # validate the problem
    validate_cse_problem(cse_problem)

    # generate the data
    u = generate_cse_data(cse_problem)

    # call the underlying implementation
    result = compute_cse(cse_problem, u)

    return result
end
