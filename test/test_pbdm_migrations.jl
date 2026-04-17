import PhysiologicallyBasedDemographicModels
using StructuredPopulationCore: StateBlockLayout, blocknames

const PBDM = PhysiologicallyBasedDemographicModels

function _population_layout(pop::PBDM.Population)
    return StateBlockLayout((stage.name => stage.delay.k for stage in pop.stages)...)
end

function _flatten_population_state(pop::PBDM.Population)
    state = Float64[]
    for stage in pop.stages
        append!(state, stage.delay.W)
    end
    return state
end

function _set_population_state!(pop::PBDM.Population, state::AbstractVector, layout::StateBlockLayout)
    length(state) == length(layout) || throw(DimensionMismatch(
        "state has length $(length(state)); expected $(length(layout))"))
    for stage in pop.stages
        range = layout[stage.name]
        stage.delay.W .= state[range]
    end
    return pop
end

function _substage_names(pop::PBDM.Population)
    names = Symbol[]
    for stage in pop.stages
        for idx in 1:stage.delay.k
            push!(names, Symbol(stage.name, "_", idx))
        end
    end
    return names
end

function _stage_totals(state::AbstractVector, layout::StateBlockLayout)
    return [sum(state[layout[name]]) for name in blocknames(layout)]
end

function _component_stage_totals(sol::CoupledMPMSolution, name::Symbol, layout::StateBlockLayout)
    totals = zeros(Float64, length(blocknames(layout)), length(sol.u))
    for (idx, snapshot) in enumerate(sol.u)
        totals[:, idx] .= _stage_totals(snapshot[name], layout)
    end
    return totals
end

function _population_step_matrix(pop_template::PBDM.Population,
        weather_day::PBDM.DailyWeather; stage_stress=nothing)
    layout = _population_layout(pop_template)
    n = length(layout)
    A = zeros(Float64, n, n)
    for col in 1:n
        pop = deepcopy(pop_template)
        basis = zeros(Float64, n)
        basis[col] = 1.0
        _set_population_state!(pop, basis, layout)
        PBDM.step_population!(pop, weather_day; stage_stress=stage_stress)
        A[:, col] .= _flatten_population_state(pop)
    end
    return A
end

function _hybrid_reference_metrics(pop_template::PBDM.Population,
        weather::PBDM.WeatherSeries,
        model::PBDM.CoupledPBDMModel,
        tspan::Tuple{Int, Int})
    pop = deepcopy(pop_template)
    ratios = Float64[]
    stresses = Vector{Vector{Float64}}()
    for day in tspan[1]:tspan[2]
        result = PBDM.step_system!(pop, PBDM.get_weather(weather, day), model)
        push!(ratios, result.supply_demand)
        push!(stresses, collect(result.stage_stress))
    end
    return (ratios = ratios, stresses = stresses)
end

@testset "PBDM pilot migrations" begin
    @testset "Exact discrete coupled migration" begin
        dev = PBDM.LinearDevelopmentRate(10.0, 35.0)

        function make_population(name, juvenile0, adult0)
            return PBDM.Population(name, [
                PBDM.LifeStage(:juvenile, PBDM.DistributedDelay(3, 45.0; W0=juvenile0), dev, 0.01),
                PBDM.LifeStage(:adult, PBDM.DistributedDelay(2, 30.0; W0=adult0), dev, 0.005),
            ])
        end

        wild = make_population(:wild, 20.0, 6.0)
        sterile = make_population(:sterile, 0.0, 0.0)

        days = [PBDM.DailyWeather(25.0, 20.0, 30.0; radiation=18.0, photoperiod=13.5)
                for _ in 1:8]
        weather = PBDM.WeatherSeries(days; day_offset=1)

        ref_rules = PBDM.AbstractInteractionRule[
            PBDM.ReproductionRule(:wild,
                (sys, w, day, p) -> 0.15 * PBDM.delay_total(sys[:wild].population.stages[2].delay))
        ]
        ref_events = PBDM.AbstractScheduledEvent[
            PBDM.SingleDayRelease(:sterile, 4.0, 3)
        ]
        ref_observables = [PBDM.Observable(:total, (sys, w, day, p) -> PBDM.total_population(sys))]

        ref_prob = PBDM.PBDMProblem(
            PBDM.MultiTypePBDM(),
            PBDM.PopulationSystem(:wild => deepcopy(wild), :sterile => deepcopy(sterile)),
            weather,
            (1, 8);
            rules = ref_rules,
            events = ref_events,
            observables = ref_observables,
        )
        ref_sol = solve(ref_prob, DirectIteration())

        wild_layout = _population_layout(wild)
        sterile_layout = _population_layout(sterile)

        mpm_prob = CoupledMPMProblem(
            PopulationSystem(
                :wild => PopulationComponent(
                    _population_step_matrix(wild, days[1]),
                    _flatten_population_state(wild);
                    stage_names = _substage_names(wild),
                    species = :fly,
                    type = :wild,
                ),
                :sterile => PopulationComponent(
                    _population_step_matrix(sterile, days[1]),
                    _flatten_population_state(sterile);
                    stage_names = _substage_names(sterile),
                    species = :fly,
                    type = :sterile,
                ),
            ),
            (1, 9);
            p = (layouts = Dict(:wild => wild_layout, :sterile => sterile_layout),),
            rules = [
                ReproductionRule(
                    :wild,
                    (sys, day, p) -> 0.15 * sum(sys[:wild].population[p.layouts[:wild][:adult]]);
                    name = :rule_1,
                ),
            ],
            events = [
                SingleDayRelease(:sterile, 4.0, 3; name = :SingleDayRelease),
            ],
            observables = [
                Observable(:total, (sys, day, p) -> total_population(sys)),
            ],
        )
        mpm_sol = solve(mpm_prob, DirectIteration())

        @test mpm_sol.retcode == :Success
        @test mpm_sol[:wild][2:end] ≈ ref_sol[:wild] atol=1e-10
        @test mpm_sol[:sterile][2:end] ≈ ref_sol[:sterile] atol=1e-10
        @test _component_stage_totals(mpm_sol, :wild, wild_layout)[:, 2:end] ≈
            ref_sol.component_stage_totals[:wild] atol=1e-10
        @test _component_stage_totals(mpm_sol, :sterile, sterile_layout)[:, 2:end] ≈
            ref_sol.component_stage_totals[:sterile] atol=1e-10
        @test mpm_sol.observables[:total][2:end] ≈ ref_sol.observables[:total] atol=1e-10
        @test mpm_sol.event_log == ref_sol.event_log
        @test getproperty.(mpm_sol.rule_log[:rule_1], :offspring) ≈
            getproperty.(ref_sol.rule_log[:rule_1], :offspring) atol=1e-10
    end

    @testset "Exact hybrid allocation migration" begin
        dev = PBDM.LinearDevelopmentRate(12.0, 35.0)
        fr = PBDM.FraserGilbertResponse(0.7)
        resp = PBDM.Q10Respiration(0.016, 2.3, 25.0)
        bdf = PBDM.BiodemographicFunctions(dev, fr, resp; label = :cotton_hybrid)
        pool = PBDM.MetabolicPool(1.0, [0.8, 1.2], [:leaf, :fruit])
        hybrid = PBDM.CoupledPBDMModel(bdf, pool; label = :cotton_hybrid)

        plant = PBDM.Population(:cotton, [
            PBDM.LifeStage(:leaf, PBDM.DistributedDelay(4, 120.0; W0=12.0), dev, 0.001),
            PBDM.LifeStage(:fruit, PBDM.DistributedDelay(3, 90.0; W0=3.0), dev, 0.002),
        ])

        weather_days = [
            PBDM.DailyWeather(24.0, 20.0, 29.0; radiation = radiation)
            for radiation in (6.0, 8.0, 5.0, 9.0, 7.0, 10.0)
        ]
        weather = PBDM.WeatherSeries(weather_days; day_offset=1)

        ref_prob = PBDM.PBDMProblem(hybrid, deepcopy(plant), weather, (1, 6))
        ref_sol = solve(ref_prob, DirectIteration())
        ref_metrics = _hybrid_reference_metrics(plant, weather, hybrid, (1, 6))

        layout = _population_layout(plant)
        stage_order = collect(blocknames(layout))
        n_stage_blocks = length(stage_order)

        mpm_prob = CoupledMPMProblem(
            PopulationSystem(
                :cotton => PopulationComponent(
                    (sys, day, p) -> _population_step_matrix(
                        p.template,
                        PBDM.get_weather(p.weather, day);
                        stage_stress = get_state(sys, :stage_stress),
                    ),
                    _flatten_population_state(plant);
                    stage_names = _substage_names(plant),
                    species = :cotton,
                );
                state = (
                    gross_supply = 0.0,
                    respiration = 0.0,
                    net_supply = 0.0,
                    demand = 0.0,
                    ratio = 1.0,
                    stage_stress = zeros(Float64, n_stage_blocks),
                    allocations = zeros(Float64, n_stage_blocks),
                ),
            ),
            (1, 7);
            p = (
                template = deepcopy(plant),
                weather = weather,
                layout = layout,
                model = hybrid,
                stage_order = stage_order,
            ),
            substeps = [
                CustomSubstep(:budget, (sys, day, p) -> begin
                    weather_day = PBDM.get_weather(p.weather, day)
                    state = sys[:cotton].population
                    stage_totals = _stage_totals(state, p.layout)
                    demand = sum(stage_totals)
                    gross_supply = PBDM.acquire(p.model.bdf.acquisition, weather_day.radiation, demand)
                    respiration = PBDM.respiration_rate(p.model.bdf.respiration, weather_day.T_mean) * demand
                    net_supply = max(0.0, gross_supply - respiration)
                    set_state!(sys, :gross_supply, gross_supply)
                    set_state!(sys, :respiration, respiration)
                    set_state!(sys, :net_supply, net_supply)
                    set_state!(sys, :demand, demand)
                    return (
                        gross_supply = gross_supply,
                        respiration = respiration,
                        net_supply = net_supply,
                        demand = demand,
                    )
                end),
                CustomSubstep(:allocation, (sys, day, p) -> begin
                    stage_totals = _stage_totals(sys[:cotton].population, p.layout)
                    demands = if length(p.model.allocation.demands) == length(stage_totals)
                        p.model.allocation.demands .* stage_totals
                    else
                        stage_totals
                    end
                    alloc_pool = PriorityAllocationPool(
                        get_state(sys, :net_supply),
                        demands,
                        p.model.allocation.labels,
                    )
                    allocations, stress = allocation_stress(alloc_pool)
                    ratio = supply_demand_index(alloc_pool)
                    set_state!(sys, :allocations, allocations)
                    set_state!(sys, :stage_stress, stress)
                    set_state!(sys, :ratio, ratio)
                    return (
                        allocations = allocations,
                        stress = stress,
                        ratio = ratio,
                    )
                end),
            ],
            observables = [
                Observable(:ratio, (sys, day, p) -> get_state(sys, :ratio)),
            ],
        )
        mpm_sol = solve(mpm_prob, DirectIteration())

        @test mpm_sol.retcode == :Success
        @test mpm_sol[:cotton][2:end] ≈ vec(sum(ref_sol.stage_totals; dims=1)) atol=1e-10
        @test _component_stage_totals(mpm_sol, :cotton, layout)[:, 2:end] ≈
            ref_sol.stage_totals atol=1e-10
        @test getproperty.(mpm_sol.substep_log[:allocation], :ratio) ≈ ref_metrics.ratios atol=1e-10
        @test mpm_sol.observables[:ratio][2:end] ≈ ref_metrics.ratios atol=1e-10
        for (actual, expected) in zip(getproperty.(mpm_sol.substep_log[:allocation], :stress), ref_metrics.stresses)
            @test actual ≈ expected atol=1e-10
        end
    end
end
