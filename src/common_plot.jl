
@userplot CSEPlot

@recipe function f(
    cp::CSEPlot,
)
    sol = cp.args[1]
    if !(sol isa SymmetricCSESolution)
        throw(ArgumentError("CSEPlot expects a SymmetricCSESolution as its first argument"))
    end
    # TODO: expand to work with asym too

    # options
    cse_label = get(plotattributes, :cse_label, "CSE")
    add_knots = get(plotattributes, :add_knots, true)
    knot_label = get(plotattributes, :knot_label, "Knots")
    add_bne = get(plotattributes, :add_bne, true)
    bne_label = get(plotattributes, :bne_label, "BNE")

    # default x and y axis labels
    xlabel --> "Private-Values"
    ylabel --> "Bids"

    # default title
    c1_string = isnan(sol.c_1) ? "" : @sprintf("; C_1=%.2e", sol.c_1)
    ftr = @sprintf("MSE=%.2e%s; C_2=%.2e", sol.mse, c1_string, sol.c_2)
    titlefontsize --> 10
    diststring = simplify_distribution_string(repr(sol.problem.distribution))
    title --> "CSE and BNE for n=$(sol.n), np=$(sol.problem.np), mc=$(sol.problem.mc)\n$(diststring)\n$(ftr)"

    # plot the CSE
    @series begin
        seriescolor := 1
        label := cse_label
        sol.cse."x", sol.cse."CSE(x)"
    end

    # optionally plot the knot points
    if (add_knots === true) && (length(sol.knot."knot(l)") > 0)
        @series begin
            seriescolor := 1
            seriestype := :scatter
            label := knot_label
            [sol.knot."knot(l-1)"[1]; sol.knot."knot(l)"], [0.0; sol.knot."CSE[knot(l)]"]
        end
    end

    # optionally plot the BNE
    if add_bne
        @series begin
            seriescolor := 2
            label := bne_label
            sol.cse."x", sol.cse."BNE(x)"
        end
    end
end
