@testset verbose = true "Test symmetric compute with 2 n values" begin
    cse_prob = SymmetricCSEProblem(inin=2, maxn=4)
    validate_cse_problem(cse_prob)
    solutions = compute_cse(cse_prob)
    @test length(solutions) == 2

    sol = solutions[1]
    @test sol.success
    @test isnan(sol.c_1)
    @test isapprox(sol.c_2, 3.715063021296306e-14)
    @test isapprox(sol.mse, 0.00020157222130046955)

    sol = solutions[2]
    @test sol.success
    @test isapprox(sol.c_1, 0.0009946887825929547)
    @test isapprox(sol.c_2, 3.87197277083241e-17)
    @test isapprox(sol.mse, 5.4743805339204674e-5)
end

@testset verbose = true "Test validate symmetric problem" begin
    @testset verbose = true "Test validate default problem" begin
        prob = SymmetricCSEProblem()
        @test validate_cse_problem(prob) === nothing
    end

    @testset verbose = true "Test validate bad distribution" begin
        prob = SymmetricCSEProblem(distribution = Normal())
        @test_throws "Only Beta distributions are supported currently" validate_cse_problem(prob)
    end
end
