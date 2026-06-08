module MatrixProjectionModels

using CommonSolve
using Distributions
using LinearAlgebra
using StructuredPopulationCore
using Random
using RecipesBase
using SciMLBase
using Statistics

# --- Layer 1: Utility functions ---
include("utils.jl")
# area_under_curve re-exported from StructuredPopulationCore

# --- Layer 2: Core types ---
include("types.jl")
export StageClassType, ActiveStage, PropaguleStage, DormantStage
export StageClass, MatrixProjectionModel, n_stages
export AbstractMPMStructure, LeslieMPM, LefkovitchMPM
# Re-export shared types from StructuredPopulationCore
export AbstractProjectionStructure
export AbstractDensityDependence, DensityIndependent, DensityDependent
export AbstractStochasticity, Deterministic, StochasticKernelResampled, StochasticParameterResampled
export DirectIteration, EigenAnalysis

# --- Layer 3: Problem/Solution ---
include("problems.jl")
export MPMProblem, remake

include("solve.jl")
using CommonSolve: solve
export solve, MPMSolution

# --- Layer 4: Coupled/stateful discrete systems ---
include("coupled.jl")
export PopulationComponent, PopulationSystem
export get_state, set_state!, has_state
export total_population, component_total, component_totals
export by_species, by_type, by_patch
export inject!, remove_fraction!
export PriorityAllocationPool, allocate, allocation_stress, supply_demand_index
export AbstractDailySubstep, apply_substep!
export StateUpdateSubstep, CustomSubstep
export StageFlux, StageTransfer, RuleEffect
export AbstractTransitionRule, AbstractScheduledEvent
export apply_rule, apply_event!
export TransferRule, ReproductionRule, MortalityRule, CustomRule
export PulseRelease, SingleDayRelease, SprayEvent, CustomEvent
export Observable
export CoupledMPMProblem, CoupledMPMSolution
export StateDependentMPMProblem, StateDependentMPMSolution
export HybridMPMProblem, HybridMPMSolution

# --- Layer 5: Analysis ---
include("analysis.jl")
# Re-export shared analysis from StructuredPopulationCore
export AbstractProjectionSolution
export eigenanalysis_power, eigenanalysis_full
export lambda, stable_distribution, reproductive_value
export sensitivity, elasticity, damping_ratio
export stochastic_growth_rate, mean_kernel, mean_matrix
export area_under_curve
# Re-export LTRE analysis
export LTREResult, StochasticLTREResult, SNALTREResult
export ExactLTREResult, RandomLTREResult
export ltre, stochastic_ltre, sna_ltre
export ltre_random, exact_ltre, gmatrix

# --- Layer 6: Matrix properties ---
include("properties.jl")
export is_ergodic, is_irreducible, is_primitive, is_leslie

# --- Layer 7: Construction (mpmsim) ---
include("construction/mortality_models.jl")
export AbstractMortalityModel
export GompertzMortality, GompertzMakehamMortality, ExponentialMortality
export SilerMortality, WeibullMortality, WeibullMakehamMortality
export model_survival

include("construction/fecundity_models.jl")
export AbstractFecundityModel
export LogisticFecundity, StepFecundity, VonBertalanffyFecundity
export NormalFecundity, HadwigerFecundity
export model_fecundity

include("construction/leslie.jl")
export make_leslie_mpm

include("construction/lefkovitch.jl")
export rand_lefko_mpm

include("construction/sets.jl")
export rand_leslie_set, rand_lefko_set

include("construction/transitions.jl")

include("construction/sampling_error.jl")
export add_mpm_error

include("density_response.jl")
export AbstractDensityResponse
export RickerDensity, BevertonHoltDensity, UsherDensity
export LogisticDensity, ThetaLogisticDensity, ConstantDensity
export DensityVitalRateSpec, apply_density

include("construction/error_estimation.jl")
export calculate_errors, compute_ci

# --- Layer 8: Life tables ---
include("life_tables/conversions.jl")
export lx_to_px, lx_to_hx, px_to_lx, px_to_hx, hx_to_lx, hx_to_px

include("life_tables/age_from_stage.jl")
export mpm_to_lx, mpm_to_px, mpm_to_hx, mpm_to_mx, mpm_to_table

# --- Layer 9: Vital rates ---
include("vital_rates/extraction.jl")
export vr_vec_survival, vr_vec_growth, vr_vec_shrinkage, vr_vec_stasis
export vr_vec_reproduction, vr_vec_dorm_enter, vr_vec_dorm_exit

include("vital_rates/averaging.jl")
export vr_survival, vr_growth, vr_shrinkage, vr_stasis, vr_fecundity

# --- Layer 10: Life history traits ---
include("life_history/life_expectancy.jl")
export life_expect_mean, life_expect_var

include("life_history/maturity.jl")
export mature_prob, mature_age, mature_distrib

include("life_history/net_repro_rate.jl")
export net_repro_rate

include("life_history/generation_time.jl")
export gen_time

include("life_history/longevity.jl")
export longevity

include("life_history/entropy.jl")
export entropy_k, entropy_k_age, entropy_k_stage, entropy_d

include("life_history/shape.jl")
export shape_surv, shape_rep

# --- Layer 11: Perturbation analysis ---
include("perturbation/matrix_perturbation.jl")
export perturb_matrix

include("perturbation/vr_perturbation.jl")
export perturb_vr

include("perturbation/transition_perturbation.jl")
export perturb_trans

include("perturbation/stochastic_perturbation.jl")
export perturb_stochastic

# --- Layer 12: Transformation ---
include("transformation/stages.jl")
export repro_stages, standard_stages, name_stages

include("transformation/split.jl")
export mpm_split

include("transformation/collapse.jl")
export mpm_collapse

include("transformation/standardize.jl")
export mpm_standardize, mpm_rearrange

# --- Layer 13: Time-lagged models ---
include("time_lag.jl")
export LaggedMPM
# Re-export time-lag functions from StructuredPopulationCore
export TimeLagStructure, expand_lag_matrix, extract_lag_components
export augment_population, extract_population
export net_repro_rate_lagged, generation_time_lagged

# --- Layer 14: Historical conditional extraction ---
include("historical_conditional.jl")
export conditional_ahistorical, conditional_difference, history_dependence

# --- Layer 15: Supplement/overwrite ---
include("supplement.jl")
export TransitionSupplement, supplement!, supplement

# --- Layer 16: SciML interface ---
include("sciml_interface.jl")
export to_discrete_problem

# --- Layer 17: Plotting ---
include("plotting.jl")

# --- Layer 18: Quasi-extinction & Markov (re-export from Core) ---
export QuasiExtinctionResult, quasi_extinction
export MarkovEnvironment, sample_next, sample_initial
export simulate_environments, project_markov

# --- Layer 19: COMPADRE extension stubs ---
include("compadre_stubs.jl")
export AbstractCompadreDB
export cdb_fetch, cdb_load, cdb_save
export cdb_matA, cdb_matU, cdb_matF, cdb_matC, cdb_metadata, cdb_id
export cdb_flag, cdb_collapse, cdb_rbind, cdb_flatten, cdb_subset, cdb_build_cdb
export mpm_mean, mpm_sd, mpm_median
export mpm_has_prop, mpm_has_active, mpm_has_dorm, mpm_first_active

end # module MatrixProjectionModels
