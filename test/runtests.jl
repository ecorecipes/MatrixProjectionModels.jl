using Test
using MatrixProjectionModels
using LinearAlgebra
using Random
using Statistics

@testset "MatrixProjectionModels.jl" begin
    include("test_types.jl")
    include("test_analysis.jl")
    include("test_properties.jl")
    include("test_solve.jl")
    include("test_leslie.jl")
    include("test_lefkovitch.jl")
    include("test_mortality_models.jl")
    include("test_fecundity_models.jl")
    include("test_sampling_error.jl")
    include("test_life_tables.jl")
    include("test_vital_rates.jl")
    include("test_life_history.jl")
    include("test_perturbation.jl")
    include("test_transformation.jl")
    include("test_transitions.jl")
    include("test_time_lag.jl")
end
