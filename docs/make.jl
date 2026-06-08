using Documenter
using MatrixProjectionModels
using StructuredPopulationCore

makedocs(;
    modules = [MatrixProjectionModels, StructuredPopulationCore],
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
        "Tutorials" => [
            "Introduction to Matrix Projection Models" => "tutorials/01_introduction.md",
            "Age-Structured (Leslie) Models" => "tutorials/02_leslie_models.md",
            "Model Construction and Random Generation" => "tutorials/03_construction.md",
            "Vital Rate Extraction and Decomposition" => "tutorials/04_vital_rates.md",
            "Life Tables and Age-from-Stage Analysis" => "tutorials/05_life_tables.md",
            "Life History Traits" => "tutorials/06_life_history.md",
            "Perturbation Analysis" => "tutorials/07_perturbation.md",
            "Stochastic Matrix Projection Models" => "tutorials/08_stochastic.md",
            "Population Simulation and Density Dependence" => "tutorials/09_simulation.md",
            "Comparative Demography with COMPADRE and COMADRE" => "tutorials/10_comparative.md",
            "Sparse Transition Constructors" => "tutorials/11_transitions.md",
            "Time-Lagged Matrix Projection Models" => "tutorials/12_time_lag.md",
            "Stage Taxonomy and Type Hierarchy" => "tutorials/13_stage_taxonomy.md",
            "Time-Lag Models and Spectral Diagnostics" => "tutorials/14_lag_and_diagnostics.md",
            "Coupled Population Systems" => "tutorials/15_coupled_systems.md",
            "Substeps, Rules, and Scheduled Events" => "tutorials/16_substeps_rules_events.md",
            "Mortality, Fecundity, and Error Estimation" => "tutorials/17_mortality_fecundity.md",
            "COMPADRE/COMADRE Database Interface" => "tutorials/18_compadre.md",
            "Population Dynamics of *Diaphorina citri* Under Seasonal Forcing" => "tutorials/19_dcitri_dynamics.md",
            "Life Table Response Experiments (LTRE)" => "tutorials/20_ltre.md",
            "Quasi-Extinction Analysis" => "tutorials/21_quasi_extinction.md",
            "Markov Environment Switching" => "tutorials/22_markov_environment.md",
            "Density-Dependent Projection" => "tutorials/23_density_dependence.md",
            "Historical Conditional Matrices" => "tutorials/24_historical_conditional.md",
        ],
    ],
)

deploydocs(;
    repo = "github.com/ecorecipes/MatrixProjectionModels.jl.git",
    push_preview = true,
)
