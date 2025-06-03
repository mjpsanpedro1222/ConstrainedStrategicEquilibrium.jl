
"""
$(TYPEDSIGNATURES)

Generate data to use for computing the asymmetric CSE case.
"""
function generate_cse_data(cse_problem::AsymmetricCSEProblem)
    # arrays for computing/storing data
    u = Array{Float64,2}(undef, cse_problem.mc, cse_problem.np)
    uni = Array{Float64,2}(undef, cse_problem.mc, cse_problem.np)

    # loop and create the data
    player_means = zeros(Float64, cse_problem.np)
    for m in 1:cse_problem.mc
        # loop over the players
        for j in 1:cse_problem.np
            # generate a uniform random number between 0 and 1
            uni[m, j] = rand(cse_problem.rng)

            # inverse beta cdf for each player
            u[m, j] = quantile(cse_problem.distributions[j], uni[m, j])
            player_means[j] += u[m, j]
        end
    end

    # computing mean and std for each player
    player_means ./= cse_problem.mc
    player_stds = zeros(Float64, cse_problem.np)
    for m = 1:cse_problem.mc
        for j = 1:cse_problem.np
            player_stds[j] += (u[m, j] - player_means[j])^2
        end
    end
    player_stds ./= cse_problem.mc
    for j in 1:cse_problem.np
        @info "mean and std player $j: $(player_means[j]), $(sqrt(player_stds[j]))"
    end

    return u
end
