
@recipe function f(sol::SymmetricCSESolution)
    xlabel --> "Private-Values"
    ylabel --> "Bids"
    label := "CSE"
    c1_string = isnan(sol.c_1) ? "" : @sprintf("; C_1=%.2e", sol.c_1)
    ftr = @sprintf("MSE=%.2e%s; C_2=%.2e", sol.mse, c1_string, sol.c_2)
    titlefontsize --> 10
    title --> "CSE and BNE for n=$(sol.n), symmetric case, $(simplify_distribution_string(repr(sol.problem.distribution)))\n$(ftr)"

    @series begin
        seriescolor := 1
        label := "CSE"
        sol.cse."x", sol.cse."CSE(x)"
    end

    @series begin
        seriescolor := 1
        seriestype := :scatter
        label := "Knots"
        [sol.knot."knot(l-1)"[1]; sol.knot."knot(l)"], [0.0; sol.knot."CSE[knot(l)]"]
    end

    @series begin
        seriescolor := 2
        label := "BNE"
        sol.cse."x", sol.cse."BNE(x)"
    end
end
