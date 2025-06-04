
function Base.show(io::IO, obj::SymmetricCSEProblem)
    print(io, "SymmetricCSEProblem(np=$(obj.np), mc=$(obj.mc), n=$(obj.inin)..$(obj.maxn), $(simplify_distribution_string(repr(obj.distribution))))")
end


# TODO: either have a generic solution type or an abstract type above them...
"""
$(TYPEDEF)

Contains the solution to the CSE problem.

$(TYPEDFIELDS)
"""
@kwdef mutable struct SymmetricCSESolution
    "A `DataFrame` containing the CSE and BNE evaluated at the given points"
    cse::DataFrame = DataFrame("x" => Float64[], "CSE(x)" => Float64[], "BNE(x)" => Float64[])
    "Mean squared error of the CSE compared to the BNE"
    mse::Float64 = NaN
    "Norm of the derivatives"
    resid::Float64 = NaN
    knot::DataFrame = DataFrame("l" => Int[], "knot(l-1)" => Float64[], "knot(l)" => Float64[], "CSE[knot(l)]" => Float64[])
    alph_bet::DataFrame = DataFrame("l" => Int[], "alph(l)" => Float64[], "bet(l)" => Float64[])
    x_f::DataFrame = DataFrame("l" => Int[], "x(l)" => Float64[], "f(l)" => Float64[])
    "Whether the calculation was successful or not"
    success::Bool = false
    "Number of function evaluations"
    nfeval::Int = 0
    "Number of Jacobians created during the solve"
    njacs::Int = 0
    "Number of factorisations of the Jacobian required for the solve"
    nfactors::Int = 0
    "Number of linear solves required for the solve"
    nsolve::Int = 0
    "Total number of iterations for the nonlinear solver"
    nsteps::Int = 0
    "Stop criteria C_1 (comparison of CSE with previous n value)"
    c_1::Float64 = NaN
    "Stop criteria C_2 (norm of the residual)"
    c_2::Float64 = NaN
    "The problem that this solution was generated for"
    problem::SymmetricCSEProblem
    "The value of n used in this solution"
    n::Int
end


function Base.show(io::IO, obj::SymmetricCSESolution)
    if obj.success
        print(io, @sprintf("SymmetricCSESolution(n=%02d, MSE=%.2e, C_1=%.2e, C_2=%.2e)", obj.n, obj.mse, obj.c_1, obj.c_2))
    else
        print(io, @sprintf("SymmetricCSESolution(n=%02d, failed)", obj.n))
    end
end
