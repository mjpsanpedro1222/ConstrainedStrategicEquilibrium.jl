"""
$(TYPEDSIGNATURES)

Function to calculate the "analytical" Bayes-Nash Equilibrium.
"""
function compute_bne(v::Float64, dist::UnivariateDistribution, N::Int)
    num = first.(quadgk.(v -> cdf(dist, v) .^ (N .- 1), 0, v))
    den = cdf(dist, v) .^ (N .- 1)
    bne = v .- num ./ den

    return bne
end
