@testset verbose = true "Test afprogs symmetric reference run for n=2..4" begin
    cse_prob = SymmetricAfrprogsCSEProblem(inin=2, maxn=4, np=4)
    validate_cse_problem(cse_prob)
    solutions = compute_cse(cse_prob)
    @test length(solutions) == 2
    @test solutions[1].success
    @test isnan(solutions[1].c_1)
    @test solutions[1].c_2 < 1e-10
    @test solutions[2].mse < 1e-3
    @test solutions[2].success
    @test !isnan(solutions[2].c_1)
    @test solutions[2].c_2 < 1e-10
    @test solutions[2].c_2 < solutions[1].c_2
    @test solutions[2].mse < 1e-4

    # load the reference data
    refdata_path = joinpath(@__DIR__, "data", "afrprogs-sym-cse-reference.txt")
    lines = readlines(refdata_path)
    data = [parse.(Float64, split(replace(line, 'D' => 'E'))) for line in lines]
    array_data = reduce(vcat, [reshape(row, 1, :) for row in data])
    cse_ref = DataFrame(array_data, ["x", "CSE(x)"])

    # compare generated CSE to reference data
    cse = solutions[1].cse[:, ["x", "CSE(x)"]]
    rtol = 1e-6
    atol = 1e-2
    @test size(cse) == size(cse_ref)
    if size(cse) != size(cse_ref)
        @error "CSE dataframes have different sizes."
    end
    for col in ("x", "CSE(x)")
        is_close = isapprox.(cse[!, col], cse_ref[!, col]; atol=atol, rtol=rtol)
        @test all(is_close)
        if !all(is_close)
            bad_rows = findall(!, is_close)
            @error "CSE dataframe differs for column `$col` at rows: " bad_rows
        end
    end
end

@testset verbose = true "Test validate symmetric problem" begin
    @testset verbose = true "Test validate default problem" begin
        prob = SymmetricAfrprogsCSEProblem()
        @test validate_cse_problem(prob) === nothing
    end
end
