@testset verbose = true "Test validate afrprogs asymmetric problem" begin
    @testset verbose = true "Test validate default problem" begin
        prob = AsymmetricAfrprogsCSEProblem()
        @test validate_cse_problem(prob) === nothing
    end

    @testset verbose = true "Test validate bad distribution" begin
        prob = AsymmetricAfrprogsCSEProblem(np=2, distributions=[Normal(), Beta()])
        @test_throws "Only Beta distributions are supported currently" validate_cse_problem(prob)

        prob = AsymmetricAfrprogsCSEProblem(np=2, distributions=[Normal(), Normal()])
        @test_throws "Only Beta distributions are supported currently" validate_cse_problem(prob)

        prob = AsymmetricAfrprogsCSEProblem(np=2, distributions=[Beta(), Beta()])
        @test validate_cse_problem(prob) === nothing
    end

    @testset verbose = true "Test validate 4 player distributions" begin
        prob = AsymmetricAfrprogsCSEProblem(np=4, distributions=[Beta(), Beta(), Beta(), Beta()])
        @test validate_cse_problem(prob) === nothing

        prob = AsymmetricAfrprogsCSEProblem(np=4, distributions=[Beta(2, 2), Beta(2, 2), Beta(3, 3), Beta(3, 3)])
        @test validate_cse_problem(prob) === nothing

        prob = AsymmetricAfrprogsCSEProblem(np=4, distributions=[Beta(2, 2), Beta(2, 3), Beta(3, 3), Beta(3, 3)])
        @test_throws "Bidders 1 and 2 must have the same distribution and the same goes for bidders 3 and 4" validate_cse_problem(prob)

        prob = AsymmetricAfrprogsCSEProblem(np=4, distributions=[Beta(2, 2), Beta(2, 2), Beta(3, 3), Beta(3, 4)])
        @test_throws "Bidders 1 and 2 must have the same distribution and the same goes for bidders 3 and 4" validate_cse_problem(prob)

        prob = AsymmetricAfrprogsCSEProblem(np=4, distributions=[Beta(3, 3), Beta(4, 4)])
        @test validate_cse_problem(prob) === nothing
        @test prob.distributions[1] == Beta(3, 3)
        @test prob.distributions[2] == Beta(3, 3)
        @test prob.distributions[3] == Beta(4, 4)
        @test prob.distributions[4] == Beta(4, 4)
    end
end

@testset verbose = true "Test afrprogs asymmetric reference run" begin
    # The original fortran asym code used n=16 and the following initial guess:
    nval = 16
    xguess = Vector{Float64}(undef, 2 * nval - 1)
    xguess[1] = -2.84827109173688
    xguess[2] = -2.85209052797506
    xguess[3] = -2.85688597597939
    xguess[4] = -2.86243170560486
    xguess[5] = -2.87005433886236
    xguess[6] = -2.87972857171734
    xguess[7] = -2.89195095174828
    xguess[8] = -2.90986034638539
    xguess[9] = -2.93411212816428
    xguess[10] = -2.96990679761558
    xguess[11] = -3.03212449908756
    xguess[12] = -3.14051739544213
    xguess[13] = -3.33100693581025
    xguess[14] = -3.69184084490141
    xguess[15] = -5.61683383184477
    xguess[16] = -7.41860036191254
    xguess[17] = 4.37486841314759
    xguess[18] = 3.55315379806175
    xguess[19] = 3.09860286531008
    xguess[20] = 2.74440338704673
    xguess[21] = 2.45415125926837
    xguess[22] = 2.19413505997764
    xguess[23] = 1.94448901027584
    xguess[24] = 1.71125414593309
    xguess[25] = 1.47745377147087
    xguess[26] = 1.23478257742541
    xguess[27] = 1.01951187118377
    xguess[28] = 0.846654741348483
    xguess[29] = 0.726435662568677
    xguess[30] = 0.619443991395723
    xguess[31] = 1.12697386952317

    # Create an asymmetric CSE problem with the following:
    cse_prob = AsymmetricAfrprogsCSEProblem(
        inin=nval,
        maxn=nval,
        solver_initial_guess=xguess,
        np=4,
    )

    validate_cse_problem(cse_prob)

    # Now compute the CSE for the problem we created
    solutions = compute_cse(cse_prob)
    @test length(solutions) == 1
    @test solutions[end].success

    # get the CSE dataframe
    cse = solutions[end].cse

    # load the reference data
    refdata_path = joinpath(@__DIR__, "data", "afrprogs-asym-cse-reference.txt")
    lines = readlines(refdata_path)
    data = [parse.(Float64, split(replace(line, 'D' => 'E'))) for line in lines]
    array_data = reduce(vcat, [reshape(row, 1, :) for row in data])
    cse_ref = DataFrame(array_data, ["x", "CSE(x) 1", "CSE(x) 2"])

    # compare generated CSE to reference data
    rtol = 1e-6
    atol = 1e-3
    @test size(cse) == size(cse_ref)
    if size(cse) != size(cse_ref)
        @error "CSE dataframes have different sizes."
    end
    for col in ("x", "CSE(x) 1", "CSE(x) 2")
        is_close = isapprox.(cse[!, col], cse_ref[!, col]; atol=atol, rtol=rtol)
        @test all(is_close)
        if !all(is_close)
            bad_rows = findall(!, is_close)
            @error "CSE dataframe differs for column `$col` at rows: " bad_rows
        end
    end
end
