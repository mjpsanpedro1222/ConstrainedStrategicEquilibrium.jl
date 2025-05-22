
"""
$(TYPEDSIGNATURES)

Generate data to use for computing the symmetric CSE case.
"""
function generate_cse_data(cse_problem::SymmetricCSEProblem)
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
            u[m, j] = quantile(cse_problem.distribution, uni[m, j])
            player_means[j] += u[m, j]
        end
    end

    # computing mean and std for each player (well players 1 and 2 are combined, as are 3 and 4)
    player_means ./= cse_problem.mc
    player_stds = zeros(Float64, cse_problem.np)
    for m = 1:cse_problem.mc
        for j = 1:cse_problem.np
            player_stds[j] += (u[m, j, 1] - player_means[j])^2
        end
    end
    player_stds ./= cse_problem.mc
    for j in 1:cse_problem.np
        @info "mean and std player $j: $(player_means[j]), $(sqrt(player_stds[j]))"
    end

    # combine stats for players 1 and 2 and players 3 and 4
    #    mean1 = (player_means[1] + player_means[2]) / 2.0
    #    mean2 = (player_means[3] + player_means[4]) / 2.0
    #    std1 = 0.0
    #    std2 = 0.0
    #    for m = 1:cse_problem.mc
    #        for j = 1:cse_problem.np
    #            if j < 3
    #                std1 += (u[m, j, 1] - mean1)^2
    #            else
    #                std2 += (u[m, j, 1] - mean2)^2
    #            end
    #        end
    #    end
    #    std1 /= (2.0 * cse_problem.mc)
    #    std2 /= (2.0 * cse_problem.mc)
    #    @info "1st 2 players mean and std: $(mean1), $(sqrt(std1))"
    #    @info "2nd 2 players mean and std: $(mean2), $(sqrt(std2))"

    return u
end
