
function simplify_distribution_string(s::String)
    m = match(r"\.?([A-Za-z_]+)\{[^}]+\}\((.*)\)", s)
    return m !== nothing ? "$(m.captures[1])($(m.captures[2]))" : s
end

