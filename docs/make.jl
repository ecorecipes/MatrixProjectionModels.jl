using Documenter
using MatrixProjectionModels
using ProjectionModels

makedocs(;
    modules = [MatrixProjectionModels, ProjectionModels],
    warnonly = true,
    authors = "Simon Frost",
    sitename = "MatrixProjectionModels.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://ecorecipes.github.io/MatrixProjectionModels.jl",
    ),
    pages = [
        "Home" => "index.md",
        "API Reference" => [
            "Types" => "api/types.md",
            "Problem & Solution" => "api/problems.md",
            "Analysis" => "api/analysis.md",
            "Matrix Properties" => "api/properties.md",
            "Construction" => "api/construction.md",
            "Life Tables" => "api/life_tables.md",
            "Vital Rates" => "api/vital_rates.md",
            "Life History Traits" => "api/life_history.md",
            "Perturbation Analysis" => "api/perturbation.md",
            "Transformation" => "api/transformation.md",
            "Time-Lag Models" => "api/time_lag.md",
            "SciML Interface" => "api/sciml.md",
            "COMPADRE Integration" => "api/compadre.md",
        ],
    ],
)

deploydocs(;
    repo = "github.com/ecorecipes/MatrixProjectionModels.jl.git",
)
