@testset "Coupled/stateful MPMs" begin
    @testset "PriorityAllocationPool" begin
        pool = PriorityAllocationPool(5.0, [2.0, 4.0], [:leaf, :fruit])
        allocations = allocate(pool)
        alloc2, stress = allocation_stress(pool)

        @test allocations ≈ [2.0, 3.0]
        @test alloc2 ≈ allocations
        @test stress ≈ [0.0, 0.25]
        @test supply_demand_index(pool) ≈ 5 / 6
    end

    @testset "PopulationSystem construction and metadata" begin
        A = [0.0 1.5; 0.4 0.6]
        wild = PopulationComponent(A, [10.0, 2.0];
            stage_names = [:juvenile, :adult],
            species = :fly,
            type = :wild,
            patch = :north)
        sterile = PopulationComponent(A, [1.0, 0.0];
            stage_names = [:juvenile, :adult],
            species = :fly,
            type = :sterile,
            patch = :north)

        sys = PopulationSystem(
            :wild => wild,
            :sterile => sterile;
            state = (release_budget = 5.0,))

        @test length(sys) == 2
        @test collect(keys(sys)) == [:wild, :sterile]
        @test total_population(sys) ≈ 13.0
        @test component_total(sys, :wild) ≈ 12.0
        @test component_totals(sys)[:sterile] ≈ 1.0
        @test length(by_species(sys, :fly)) == 2
        @test length(by_type(sys, :wild)) == 1
        @test length(by_patch(sys, :north)) == 2
        @test get_state(sys, :release_budget) == 5.0

        set_state!(sys, :release_budget, 7.5)
        @test get_state(sys, :release_budget) == 7.5
        @test_throws ArgumentError PopulationSystem()
    end

    @testset "Events and rules" begin
        A = [0.0 1.0; 0.5 0.3]
        sys = PopulationSystem(
            :wild => PopulationComponent(A, [4.0, 2.0]; stage_names = [:juvenile, :adult]),
            :sterile => PopulationComponent(A, [0.0, 1.0]; stage_names = [:juvenile, :adult]))

        release = PulseRelease(:sterile, 3.0, 2; stage_idx = 2, start_day = 0)
        @test apply_event!(release, sys, 0, nothing)
        @test !apply_event!(release, sys, 1, nothing)
        @test component_total(sys, :sterile) ≈ 4.0

        transfer = TransferRule(:wild, :sterile, (sys, day, p) -> 0.5;
            source_stage = 2,
            target_stage = 1,
            name = :marking)
        effect = apply_rule(transfer, sys, 0, nothing)
        @test effect.metrics.transferred ≈ 1.0

        reproduction = ReproductionRule(:wild, (sys, day, p) -> sys[:wild].population[2] * 2.0;
            name = :wild_births)
        births = apply_rule(reproduction, sys, 0, nothing)
        @test births.metrics.offspring ≈ 4.0

        mortality = MortalityRule(:wild, (sys, day, p) -> 0.25; stage = 1, name = :juvenile_loss)
        deaths = apply_rule(mortality, sys, 0, nothing)
        @test deaths.metrics.removed ≈ 1.0
    end

    @testset "CoupledMPMProblem solve" begin
        wild_model = [0.0 2.0; 0.4 0.6]
        sterile_model = [0.0 0.0; 0.2 0.7]

        sys = PopulationSystem(
            :wild => PopulationComponent(wild_model, [10.0, 2.0];
                stage_names = [:juvenile, :adult],
                species = :fly,
                type = :wild),
            :sterile => PopulationComponent(sterile_model, [0.0, 0.0];
                stage_names = [:juvenile, :adult],
                species = :fly,
                type = :sterile);
            state = (fertility_scale = 0.5,))

        events = [SingleDayRelease(:sterile, 4.0, 0; stage_idx = 2)]
        rules = [
            ReproductionRule(:wild,
                (sys, day, p) -> get_state(sys, :fertility_scale) * sys[:wild].population[2];
                name = :wild_births),
            TransferRule(:wild, :sterile,
                (sys, day, p) -> 0.1;
                source_stage = 2,
                target_stage = 2,
                name = :sterilization)
        ]
        observables = [Observable(:total, (sys, day, p) -> total_population(sys))]

        prob = CoupledMPMProblem(sys, (0, 3);
            rules = rules,
            events = events,
            observables = observables)
        sol = solve(prob, DirectIteration())

        @test sol.retcode == :Success
        @test length(sol.t) == 4
        @test length(sol[:wild]) == 4
        @test length(sol[:sterile]) == 4
        @test length(sol.observables[:total]) == 4
        @test length(sol.event_log) == 1
        @test length(sol.rule_log[:wild_births]) == 3
        @test sol[:sterile][2] > sol[:sterile][1]
        @test sol.observables[:total][end] > 0
    end

    @testset "Dynamic component matrices and SciML lowering" begin
        dynamic = PopulationComponent(
            (sys, day, p) -> [0.0 get_state(sys, :bonus); 0.5 0.4],
            [5.0, 1.0];
            stage_names = [:juvenile, :adult])
        sys = PopulationSystem(:population => dynamic; state = (bonus = 1.5,))
        prob = CoupledMPMProblem(sys, (0, 2))
        dp = to_discrete_problem(prob)

        @test dp.u0 isa PopulationSystem
        sol = solve(prob, DirectIteration())
        @test sol[:population][2] ≈ sum([0.0 1.5; 0.5 0.4] * [5.0, 1.0])
    end

    @testset "Ordered hybrid substeps" begin
        p = (
            radiation = [6.0, 2.0],
            efficiency = 1.0,
            respiration = 0.2,
            vegetative_demand = 0.5,
            reproductive_demand = 1.5,
        )

        dynamic = PopulationComponent(
            (sys, day, p) -> begin
                stress = get_state(sys, :stress)
                ratio = get_state(sys, :ratio)
                [0.0 1.5 * ratio;
                 0.4 * (1 - stress[1]) 0.6 * (1 - stress[2])]
            end,
            [5.0, 2.0];
            stage_names = [:vegetative, :reproductive])

        sys = PopulationSystem(:plant => dynamic;
            state = (
                net_supply = 0.0,
                ratio = 1.0,
                stress = [0.0, 0.0],
                allocations = [0.0, 0.0],
            ))

        substeps = [
            StateUpdateSubstep(:budget, :net_supply,
                (sys, day, p) -> max(0.0, p.radiation[day + 1] * p.efficiency -
                    p.respiration * total_population(sys))),
            CustomSubstep(:allocation, (sys, day, p) -> begin
                demands = [
                    sys[:plant].population[1] * p.vegetative_demand,
                    sys[:plant].population[2] * p.reproductive_demand,
                ]
                pool = PriorityAllocationPool(get_state(sys, :net_supply), demands,
                    [:vegetative, :reproductive])
                allocations, stress = allocation_stress(pool)
                ratio = supply_demand_index(pool)
                set_state!(sys, :allocations, allocations)
                set_state!(sys, :stress, stress)
                set_state!(sys, :ratio, ratio)
                return (allocations = allocations, stress = stress, ratio = ratio)
            end),
        ]

        prob = CoupledMPMProblem(sys, (0, 2); p = p, substeps = substeps,
            observables = [Observable(:ratio, (sys, day, p) -> get_state(sys, :ratio))])
        sol = solve(prob, DirectIteration())

        @test sol.retcode == :Success
        @test length(sol.substep_log[:budget]) == 2
        @test length(sol.substep_log[:allocation]) == 2
        @test sol.substep_log[:budget][1].net_supply ≈ 6.0 - 0.2 * 7.0
        @test sol.substep_log[:allocation][1].ratio ≈ supply_demand_index(
            PriorityAllocationPool(6.0 - 0.2 * 7.0, [2.5, 3.0], [:vegetative, :reproductive]))
        @test sol.component_matrices[:plant][1] ≈ [
            0.0 1.5 * sol.substep_log[:allocation][1].ratio;
            0.4 * (1 - sol.substep_log[:allocation][1].stress[1]) 0.6 * (1 - sol.substep_log[:allocation][1].stress[2])
        ]
        @test sol.observables[:ratio][2] ≈ sol.substep_log[:allocation][1].ratio
    end

    @testset "Mass-balance guard" begin
        sys = PopulationSystem(
            :wild => PopulationComponent([0.0 1.0; 0.5 0.3], [1.0, 0.0];
                stage_names = [:juvenile, :adult]),
            :sink => PopulationComponent([0.0 0.0; 0.0 1.0], [0.0, 0.0];
                stage_names = [:juvenile, :adult]))

        bad_rule = CustomRule(:too_much, (sys, day, p) ->
            StageTransfer(:wild, 1, :sink, 1, 10.0))
        prob = CoupledMPMProblem(sys, (0, 1); rules = [bad_rule])

        @test_throws ArgumentError solve(prob, DirectIteration())
    end
end
