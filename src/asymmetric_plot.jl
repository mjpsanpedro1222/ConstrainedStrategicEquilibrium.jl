
@recipe function f(sol::AsymmetricCSESolution)
    xlabel --> "Private-Values"
    ylabel --> "Bids"
    label := "CSE 1"

    c1_string = ""
    if !isnan(sol.c_1[:bidder1])
        c1_string = @sprintf("(%.2e, %.2e)", sol.c_1[:bidder1], sol.c_1[:bidder2])
    end
    ftr = @sprintf("%s; C_2=%.2e", c1_string, sol.c_2)
    titlefontsize --> 10
    dist2 = length(sol.problem.distributions) == 2 ? sol.problem.distributions[2] : sol.problem.distributions[3]
    diststring = @sprintf("(1) %s; (2) %s", simplify_distribution_string(repr(sol.problem.distributions[1])), simplify_distribution_string(repr(dist2)))
    title --> "CSE for n=$(sol.n), np=$(sol.problem.np), mc=$(sol.problem.mc), asymmetric case\n$(diststring)\n$(ftr)"

    @series begin
        label := "CSE 2"
        sol.cse."x", sol.cse."CSE(x) 2"
    end

    @series begin
        seriestype := :scatter
        label := "Knots 1"
        [sol.knot[:bidder1]."knot(l-1)"[1]; sol.knot[:bidder1]."knot(l)"], [0.0; sol.knot[:bidder1]."CSE[knot(l)]"]
    end

    @series begin
        seriestype := :scatter
        label := "Knots 2"
        [sol.knot[:bidder2]."knot(l-1)"[1]; sol.knot[:bidder2]."knot(l)"], [0.0; sol.knot[:bidder2]."CSE[knot(l)]"]
    end

    # return data
    sol.cse."x", sol.cse."CSE(x) 1"
end

