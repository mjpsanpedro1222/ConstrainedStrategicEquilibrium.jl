push!(LOAD_PATH, joinpath(@__DIR__, ".."))

using Documenter
using DocumenterCitations
using Literate
using ConstrainedStrategicEquilibrium

# generate examples, based on https://github.com/CliMA/Oceananigans.jl/blob/1c2a6f8752b6425bf30d856f8ba0aa681c0ab818/docs/make.jl

const EXAMPLES_INPUT_DIR = joinpath(@__DIR__, "..", "examples")
const EXAMPLES_OUTPUT_DIR = joinpath(@__DIR__, "src", "generated")
const EXAMPLES_RELATIVE_PATH = "generated"

examples = Dict(
    "afr-progs-sym" => "Symmetric CSE - piecewise linear",
    "jae-poly-1" => "Symmetric CSE - polynomial",
)

@info "Generating examples"
example_pages = Array{Pair{String,String},1}()
for (example_name, example_title) in examples
    @info "Generating exmaple: $example_name"

    example_input_file = joinpath(EXAMPLES_INPUT_DIR, example_name * ".jl")
    example_output_file = joinpath(EXAMPLES_OUTPUT_DIR, example_name * ".md")

    # check if we need to generate the markdown
    generate_markdown = true
    if isfile(example_output_file)
        mtime_generated = mtime(example_output_file)
        mtime_source = mtime(example_input_file)

        # don't need to generate if the generated file was modified more recently than the source
        if mtime_generated > mtime_source
            @info "Skipping generation of $example_name (already generated)"
            generate_markdown = false
        end
    end

    if generate_markdown
        Literate.markdown(
            example_input_file,
            EXAMPLES_OUTPUT_DIR;
            flavor=Literate.DocumenterFlavor(),
            execute=true
        )
    end

    push!(example_pages, example_title => joinpath(EXAMPLES_RELATIVE_PATH, example_name * ".md"))
end
@info "Finished generating examples"

# organise pages

basics_pages = [
    "CSE problems" => "basics/cse_problem.md",
    "CSE solutions" => "basics/cse_solution.md",
    "Solving a CSE" => "basics/solve.md",
    "Plotting a CSE solution" => "basics/plotting.md",
]

pages = [
    "Home" => "index.md",
    "Basics" => basics_pages,
    "Examples" => example_pages,
    "API" => "api.md"
]

# build and deploy docs

format = Documenter.HTML(
    prettyurls=get(ENV, "CI", nothing) == "true",
    assets=String["assets/citations.css"],
)

# citations
bib = CitationBibliography(
    joinpath(@__DIR__, "src", "refs.bib");
    style=:numeric,
)

makedocs(
    sitename="ConstrainedStrategicEquilibrium.jl",
    format=format,
    pages=pages,
    plugins=[bib],
)

#deploydocs(
#    repo = "",
#)
