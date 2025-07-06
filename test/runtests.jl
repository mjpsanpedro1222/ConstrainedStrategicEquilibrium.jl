using ConstrainedStrategicEquilibrium
using DataFrames
using Distributions
using Test

@testset verbose = true "ConstrainedStrategicEquilibrium tests" begin
    @testset verbose = true "Symmetric CSE tests" begin
        include("test_symmetric.jl")
    end

    @testset verbose = true "Asymmetric CSE tests" begin
        include("test_afrprogs_asymmetric.jl")
    end
end
