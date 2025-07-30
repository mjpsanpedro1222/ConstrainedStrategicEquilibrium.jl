
@userplot CSEPlot

"""
    cseplot(sol::CSESolution; kwargs...)

Plot the given `CSESolution`. Some additional keyword arguments are available (in addition to the
standard plotting options, such as `dpi`, `title`, etc.).

When `sol` is a `SymmetricCSESolution` these are (default values shown):

- `cse_label::String = "CSE"` - the label for the CSE line
- `cse_colour::Symbol = :1` - the colour for the CSE line
- `add_knots::Bool = true` - whether or not to show the knot points on the graph
- `knots_label::String = "Knots"` - the label for the knot points
- `add_bne::Bool = true` - whether or not to show the BNE on the graph
- `bne_label::String = "BNE"` - the label for the BNE line
- `bne_colour::Symbol = :2` - the colour for the BNE line

When `sol` is an `AsymmetricCSESolution` the extra options are (default values shown):

- `cse_label::Dict{Symbol, String} = Dict(:bidder1 => "CSE (1), :bidder2 => "CSE (2)")` - the labels
- `cse_colour::Dict{Symbol, Symbol} = Dict(:bidder1 => :1, :bidder2 => :2` - the colours
  for the CSE lines for each bidder (the keys in the dictionary must be as shown)
- `add_knots::Bool = true` - whether or not to show the knot points on the graph
- `knots_label::Dict{Symbol, String} = Dict(:bidder1 => "Knots (1)", :bidder2 => "Knots (2)")` - the
  labels for the knot points for each bidder (the keys in the dictionary must be as shown)

In addition to the above options, `cseplot` will set a default `title` and `xlabel` and `ylabel`,
however these can be overriden with keyword arguments.
"""
function cseplot end

@recipe function f(
    cp::CSEPlot,
)
    sol = cp.args[1]
    if !(sol isa CSESolution)
        throw(ArgumentError("CSEPlot expects a CSESolution as its first argument"))
    end

    # default x and y axis labels
    xlabel --> "Private-Values"
    ylabel --> "Bids"

    # common options
    add_knots = get(plotattributes, :add_knots, true)

    if isa(sol, SymmetricCSESolution)
        # options
        cse_label = get(plotattributes, :cse_label, "CSE")
        cse_colour = get(plotattributes, :cse_colour, :1)
        knot_label = get(plotattributes, :knot_label, "Knots")
        add_bne = get(plotattributes, :add_bne, true)
        bne_label = get(plotattributes, :bne_label, "BNE")
        bne_colour = get(plotattributes, :bne_colour, :2)

        # default title
        c1_string = isnan(sol.c_1) ? "" : @sprintf("; C_1=%.2e", sol.c_1)
        ftr = @sprintf("MSE=%.2e%s; C_2=%.2e", sol.mse, c1_string, sol.c_2)
        titlefontsize --> 10
        diststring = simplify_distribution_string(repr(sol.problem.distribution))
        title --> "CSE and BNE for n=$(sol.n), np=$(sol.problem.np), mc=$(sol.problem.mc)\n$(diststring)\n$(ftr)"

        # plot the CSE
        @series begin
            seriescolor := cse_colour
            label := cse_label
            sol.cse."x", sol.cse."CSE(x)"
        end

        # optionally plot the knot points
        if (add_knots === true) && (length(sol.knot."knot(l)") > 0)
            @series begin
                seriescolor := cse_colour
                seriestype := :scatter
                label := knot_label
                [sol.knot."knot(l-1)"[1]; sol.knot."knot(l)"], [0.0; sol.knot."CSE[knot(l)]"]
            end
        end

        # optionally plot the BNE
        if add_bne
            @series begin
                seriescolor := bne_colour
                label := bne_label
                sol.cse."x", sol.cse."BNE(x)"
            end
        end
    elseif isa(sol, AsymmetricCSESolution)
        # options
        cse_label = get(plotattributes, :cse_label, Dict(:bidder1 => "CSE (1)", :bidder2 => "CSE (2)"))
        cse_colour = get(plotattributes, :cse_colour, Dict(:bidder1 => :1, :bidder2 => :2))
        knot_label = get(plotattributes, :knot_label, Dict(:bidder1 => "Knots (1)", :bidder2 => "Knots (2)"))

        # default title
        c1_string = ""
        if !isnan(sol.c_1[:bidder1])
            c1_string = @sprintf("C_1=(%.2e, %.2e); ", sol.c_1[:bidder1], sol.c_1[:bidder2])
        end
        ftr = @sprintf("%sC_2=%.2e", c1_string, sol.c_2)
        titlefontsize --> 10
        dist2 = length(sol.problem.distributions) == 2 ? sol.problem.distributions[2] : sol.problem.distributions[3]
        diststring = @sprintf("(1) %s; (2) %s", simplify_distribution_string(repr(sol.problem.distributions[1])), simplify_distribution_string(repr(dist2)))
        title --> "CSE for n=$(sol.n), np=$(sol.problem.np), mc=$(sol.problem.mc)\n$(diststring)\n$(ftr)"

        # plot the CSE for bidder 1
        @series begin
            seriescolor := cse_colour[:bidder1]
            label := cse_label[:bidder1]
            sol.cse."x", sol.cse."CSE(x) 1"
        end

        # optionally, plot the knots for bidder 1
        if (add_knots === true) && (length(sol.knot[:bidder1]."knot(l)") > 0)
            @series begin
                seriescolor := cse_colour[:bidder1]
                seriestype := :scatter
                label := knot_label[:bidder1]
                [sol.knot[:bidder1]."knot(l-1)"[1]; sol.knot[:bidder1]."knot(l)"], [0.0; sol.knot[:bidder1]."CSE[knot(l)]"]
            end
        end

        # plot the CSE for bidder 2
        @series begin
            seriescolor := cse_colour[:bidder2]
            label := cse_label[:bidder2]
            sol.cse."x", sol.cse."CSE(x) 2"
        end

        # optionally, plot the knots for bidder 2
        if (add_knots === true) && (length(sol.knot[:bidder2]."knot(l)") > 0)
            @series begin
                seriescolor := cse_colour[:bidder2]
                seriestype := :scatter
                label := knot_label[:bidder2]
                [sol.knot[:bidder2]."knot(l-1)"[1]; sol.knot[:bidder2]."knot(l)"], [0.0; sol.knot[:bidder2]."CSE[knot(l)]"]
            end
        end
    else
        @error "CSEPlot not implemented for $(typeof(sol))"
    end
end
