using ConstrainedStrategicEquilibrium
using Distributions
using Test

@testset verbose = true "ConstrainedStrategicEquilibrium tests" begin
    @testset verbose = true "Symmetric CSE tests" begin
        include("test_symmetric.jl")
    end
end
