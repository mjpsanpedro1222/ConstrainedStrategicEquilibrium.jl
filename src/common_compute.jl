
"""
$(TYPEDSIGNATURES)

Solve the given CSE problem.

This is a wrapper function that will first validate the problem and generate the data
before solving the CSE problem and returning the vector of solutions.

# Examples
```jldoctest
julia> prob = SymmetricJaePoly1CSEProblem()
SymmetricJaePoly1CSEProblem(np=4, mc=10000, n=1..12, Distributions.Kumaraswamy{Float64}(a=2.5, b=3.5))

julia> solutions = compute_cse(prob)
[ Info: mean and std player 1: 0.4957361972625512, 0.18748491448182267
[ Info: mean and std player 2: 0.5000858627485552, 0.18774978658869976
[ Info: mean and std player 3: 0.496651263173003, 0.18645006296589822
[ Info: mean and std player 4: 0.5021607262564288, 0.18574470364076182
[ Info: Computing: SymmetricJaePoly1CSEProblem(np=4, mc=10000, n=1..12, Kumaraswamy(a=2.5, b=3.5))
[ Info: SymmetricCSESolution(n=01, MSE=1.63e-03, C_1=NaN, C_2=9.48e-16)
[ Info: SymmetricCSESolution(n=02, MSE=6.62e-04, C_1=5.80e-03, C_2=1.91e-16)
[ Info: SymmetricCSESolution(n=03, MSE=9.80e-05, C_1=3.45e-03, C_2=1.15e-13)
[ Info: SymmetricCSESolution(n=04, MSE=1.75e-05, C_1=5.70e-04, C_2=4.09e-14)
[ Info: SymmetricCSESolution(n=05, MSE=7.30e-05, C_1=1.24e-03, C_2=1.00e-16)
[ Info: SymmetricCSESolution(n=06, MSE=6.51e-05, C_1=1.64e-03, C_2=1.70e-14)
[ Info: SymmetricCSESolution(n=07, MSE=1.95e-05, C_1=1.23e-03, C_2=4.07e-16)
[ Info: SymmetricCSESolution(n=08, MSE=1.94e-06, C_1=5.75e-04, C_2=2.06e-14)
[ Info: SymmetricCSESolution(n=09, MSE=2.59e-07, C_1=8.99e-05, C_2=2.10e-13)
[ Info: SymmetricCSESolution(n=10, MSE=1.17e-06, C_1=1.53e-04, C_2=6.20e-14)
[ Info: SymmetricCSESolution(n=11, MSE=1.53e-06, C_1=1.58e-05, C_2=3.61e-13)
[ Info: SymmetricCSESolution(n=12, MSE=5.77e-07, C_1=4.83e-05, C_2=8.63e-14)
12-element Vector{ConstrainedStrategicEquilibrium.SymmetricCSESolution}:
 SymmetricCSESolution(n=01, MSE=1.63e-03, C_1=NaN, C_2=9.48e-16)
 SymmetricCSESolution(n=02, MSE=6.62e-04, C_1=5.80e-03, C_2=1.91e-16)
 SymmetricCSESolution(n=03, MSE=9.80e-05, C_1=3.45e-03, C_2=1.15e-13)
 SymmetricCSESolution(n=04, MSE=1.75e-05, C_1=5.70e-04, C_2=4.09e-14)
 SymmetricCSESolution(n=05, MSE=7.30e-05, C_1=1.24e-03, C_2=1.00e-16)
 SymmetricCSESolution(n=06, MSE=6.51e-05, C_1=1.64e-03, C_2=1.70e-14)
 SymmetricCSESolution(n=07, MSE=1.95e-05, C_1=1.23e-03, C_2=4.07e-16)
 SymmetricCSESolution(n=08, MSE=1.94e-06, C_1=5.75e-04, C_2=2.06e-14)
 SymmetricCSESolution(n=09, MSE=2.59e-07, C_1=8.99e-05, C_2=2.10e-13)
 SymmetricCSESolution(n=10, MSE=1.17e-06, C_1=1.53e-04, C_2=6.20e-14)
 SymmetricCSESolution(n=11, MSE=1.53e-06, C_1=1.58e-05, C_2=3.61e-13)
 SymmetricCSESolution(n=12, MSE=5.77e-07, C_1=4.83e-05, C_2=8.63e-14)
```
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
