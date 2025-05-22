module ConstrainedStrategicEquilibrium

using Random
using Distributions
using NonlinearSolve
using SpecialFunctions
using LinearAlgebra
using QuadGK
using DocStringExtensions
using DataFrames
using SciMLBase
using RecipesBase
using Printf
using Logging

include("utils.jl")

include("base_types.jl")
include("bne.jl")

# include general symmetric case things
include("symmetric_data.jl")
include("symmetric_types.jl")
include("symmetric_compute.jl")
include("symmetric_plot.jl")

# symmetric afr progs case
include("symmetric_afrprogs.jl")

# symmetric jae_poly_1 case
include("symmetric_jae_poly_1.jl")

# TODO: include the asymmetric case


# define exports
export SymmetricAfrprogsCSEProblem, validate_cse_problem, compute_cse, generate_cse_data
export SymmetricJaePoly1CSEProblem

end # module ConstrainedStrategicEquilibrium
