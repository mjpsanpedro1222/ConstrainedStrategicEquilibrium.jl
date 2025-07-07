
@recipe function f(sol::AsymmetricCSESolution)
    xlabel --> "Private-Values"
    ylabel --> "Bids"
    label := "CSE 1"

    @series begin
        label := "CSE 2"
        sol.cse."x", sol.cse."CSE(x) 2"
    end

    # return data
    sol.cse."x", sol.cse."CSE(x) 1"
end

