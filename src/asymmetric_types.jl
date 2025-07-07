

function Base.show(io::IO, obj::AsymmetricCSEProblem)
    print(io, "AsymmetricCSEProblem(np=$(obj.np), mc=$(obj.mc), n=$(obj.inin)..$(obj.maxn))")
end


# TODO: either have a generic solution type or an abstract type above them...
"""
$(TYPEDEF)

Contains the solution to the asymmetric CSE problem.

$(TYPEDFIELDS)
"""
@kwdef mutable struct AsymmetricCSESolution <: CSESolution
    "A `DataFrame` containing the CSE and BNE evaluated at the given points"
    cse::DataFrame = DataFrame("x" => Float64[], "CSE(x) 1" => Float64[], "CSE(x) 2" => Float64[])
    "Norm of the derivatives"
    resid::Float64 = NaN
    "The knots for both bidders"
    knot::Dict{Symbol,DataFrame} = Dict(
        :bidder1 => DataFrame("l" => Int[], "knot(l-1)" => Float64[], "knot(l)" => Float64[], "CSE[knot(l)]" => Float64[]),
        :bidder2 => DataFrame("l" => Int[], "knot(l-1)" => Float64[], "knot(l)" => Float64[], "CSE[knot(l)]" => Float64[]),
    )
    "Alpha and beta values for the piecewise linear functions for both bidders"
    alph_bet::Dict{Symbol,DataFrame} = Dict(
        :bidder1 => DataFrame("l" => Int[], "alph(l)" => Float64[], "bet(l)" => Float64[]),
        :bidder2 => DataFrame("l" => Int[], "alph(l)" => Float64[], "bet(l)" => Float64[]),
    )
    "x values and derivatives at the final point"
    x_f::Dict{Symbol,DataFrame} = Dict(
        :bidder1 => DataFrame("l" => Int[], "x(l)" => Float64[], "f(l)" => Float64[]),
        :bidder2 => DataFrame("l" => Int[], "x(l)" => Float64[], "f(l)" => Float64[]),
    )
    "Norm of the derivatives for each player"
    norm_derivatives::Dict{Symbol,Float64} = Dict(
        :bidder1 => NaN,
        :bidder2 => NaN,
    )
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
    "Stop criteria C_1 (comparison of CSE with previous n value) for each player"
    c_1::Dict{Symbol,Float64} = Dict(
        :bidder1 => NaN,
        :bidder2 => NaN,
    )
    "Stop criteria C_2 (norm of the residual)"
    c_2::Float64 = NaN
    "The problem that this solution was generated for"
    problem::AsymmetricCSEProblem
    "The value of n used in this solution"
    n::Int
    "The data used in generating the solution"
    u::Array{Float64}
    "The solution object returned by the solver"
    solver_solution::Union{SciMLBase.NonlinearSolution,Nothing} = nothing
end


function Base.show(io::IO, obj::AsymmetricCSESolution)
    if obj.success
        print(io, @sprintf("AsymmetricCSESolution(n=%02d, C_1=(%.2e, %.2e), C_2=%.2e)", obj.n, obj.c_1[:bidder1], obj.c_1[:bidder2], obj.c_2))
    else
        print(io, @sprintf("AsymmetricCSESolution(n=%02d, failed)", obj.n))
    end
end
