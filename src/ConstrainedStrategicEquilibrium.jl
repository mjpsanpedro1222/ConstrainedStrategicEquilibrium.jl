module ConstrainedStrategicEquilibrium

using Random
using Distributions
using NonlinearSolve
import NonlinearSolve.NonlinearSolveBase: AbstractNonlinearAlgorithm
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
include("common_compute.jl")
include("common_solve.jl")

# include general symmetric case things
include("symmetric_data.jl")
include("symmetric_types.jl")
include("symmetric_plot.jl")

# symmetric afr progs case
include("symmetric_afrprogs.jl")

# symmetric jae_poly_1 case
include("symmetric_jae_poly_1.jl")

# include the asymmetric case
include("asymmetric_types.jl")
include("asymmetric_data.jl")
include("asymmetric_afrprogs.jl")
include("asymmetric_plot.jl")

# define exports
export SymmetricAfrprogsCSEProblem, validate_cse_problem, compute_cse, generate_cse_data
export SymmetricJaePoly1CSEProblem
export AsymmetricAfrprogsCSEProblem

end # module ConstrainedStrategicEquilibrium
