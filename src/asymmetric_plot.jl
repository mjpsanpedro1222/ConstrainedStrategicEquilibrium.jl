
@recipe function f(sol::AsymmetricCSESolution)
    xlabel --> "Private-Values"
    ylabel --> "Bids"
    label := "CSE 1"

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

