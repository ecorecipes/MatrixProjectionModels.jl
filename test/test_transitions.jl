@testset "Sparse transition constructors" begin

    @testset "A-only construction" begin
        mpm = MatrixProjectionModel([:seed, :small, :large],
            (:seed => :small) => 0.2,
            (:small => :large) => 0.4,
            (:small => :small) => 0.3,
            (:large => :large) => 0.7,
            (:large => :seed) => 5.0)

        @test size(mpm) == (3, 3)
        @test mpm.stage_names == [:seed, :small, :large]
        # (from => to) => val  →  A[to, from]
        @test mpm.A[2, 1] ≈ 0.2  # seed → small
        @test mpm.A[3, 2] ≈ 0.4  # small → large
        @test mpm.A[2, 2] ≈ 0.3  # small → small
        @test mpm.A[3, 3] ≈ 0.7  # large → large
        @test mpm.A[1, 3] ≈ 5.0  # large → seed
        # U = A, F = 0, C = 0
        @test mpm.U ≈ mpm.A
        @test all(mpm.F .== 0)
        @test all(mpm.C .== 0)
    end

    @testset "U/F/C decomposition" begin
        mpm = MatrixProjectionModel([:seed, :small, :large];
            U = [(:seed => :small) => 0.2, (:small => :large) => 0.4,
                 (:small => :small) => 0.3, (:large => :large) => 0.7],
            F = [(:large => :seed) => 5.0, (:small => :seed) => 1.0])

        @test size(mpm) == (3, 3)
        @test mpm.stage_names == [:seed, :small, :large]
        # U matrix
        @test mpm.U[2, 1] ≈ 0.2
        @test mpm.U[3, 2] ≈ 0.4
        @test mpm.U[2, 2] ≈ 0.3
        @test mpm.U[3, 3] ≈ 0.7
        @test mpm.U[1, 3] ≈ 0.0
        # F matrix
        @test mpm.F[1, 3] ≈ 5.0
        @test mpm.F[1, 2] ≈ 1.0
        @test all(mpm.F[2:3, :] .== 0)
        # C = 0
        @test all(mpm.C .== 0)
        # A = U + F + C
        @test mpm.A ≈ mpm.U .+ mpm.F .+ mpm.C
    end

    @testset "Stage names preserved" begin
        names = [:juvenile, :subadult, :adult]
        mpm = MatrixProjectionModel(names,
            (:juvenile => :subadult) => 0.5,
            (:adult => :juvenile) => 3.0)
        @test mpm.stage_names == names
        @test n_stages(mpm) == 3
    end

    @testset "Multiple entries to same cell sum" begin
        mpm = MatrixProjectionModel([:a, :b],
            (:a => :b) => 0.3,
            (:a => :b) => 0.2)
        @test mpm.A[2, 1] ≈ 0.5  # 0.3 + 0.2
    end

    @testset "Unknown stage name throws" begin
        @test_throws ArgumentError MatrixProjectionModel([:a, :b],
            (:a => :c) => 0.5)
        @test_throws ArgumentError MatrixProjectionModel([:a, :b],
            (:c => :a) => 0.5)
        @test_throws ArgumentError MatrixProjectionModel([:a, :b];
            U = [(:a => :c) => 0.5])
    end

    @testset "Empty transitions produce zero matrix" begin
        mpm = MatrixProjectionModel([:a, :b])
        @test all(mpm.A .== 0)
        @test size(mpm) == (2, 2)

        mpm2 = MatrixProjectionModel([:a, :b]; U=[], F=[], C=[])
        @test all(mpm2.A .== 0)
    end

    @testset "Type promotion" begin
        # Int values should be promoted to Float64
        mpm = MatrixProjectionModel([:a, :b],
            (:a => :b) => 1,
            (:b => :a) => 2)
        @test eltype(mpm) == Float64
        @test mpm.A[2, 1] ≈ 1.0
        @test mpm.A[1, 2] ≈ 2.0

        # Mixed Int and Float64
        mpm2 = MatrixProjectionModel([:a, :b];
            U = [(:a => :b) => 1],
            F = [(:b => :a) => 2.5])
        @test eltype(mpm2) == Float64
    end

    @testset "Equivalence with manual matrix construction" begin
        # Build manually
        A = [0.0 0.0 5.0;
             0.2 0.3 0.0;
             0.0 0.4 0.7]
        mpm_manual = MatrixProjectionModel(A;
            stage_names=[:seed, :small, :large])

        # Build via sparse transitions
        mpm_sparse = MatrixProjectionModel([:seed, :small, :large],
            (:seed => :small) => 0.2,
            (:small => :large) => 0.4,
            (:small => :small) => 0.3,
            (:large => :large) => 0.7,
            (:large => :seed) => 5.0)

        @test mpm_sparse.A ≈ mpm_manual.A
        @test lambda(mpm_sparse) ≈ lambda(mpm_manual)
    end

    @testset "Leslie-like model via sparse entries" begin
        survival = [0.8, 0.9]
        fecundity = [0.0, 1.5, 3.0]
        mpm_leslie = make_leslie_mpm(survival, fecundity)

        mpm_sparse = MatrixProjectionModel([:age_1, :age_2, :age_3];
            U = [(:age_1 => :age_2) => 0.8,
                 (:age_2 => :age_3) => 0.9],
            F = [(:age_1 => :age_1) => 0.0,
                 (:age_2 => :age_1) => 1.5,
                 (:age_3 => :age_1) => 3.0])

        @test mpm_sparse.A ≈ mpm_leslie.A
        @test mpm_sparse.U ≈ mpm_leslie.U
        @test mpm_sparse.F ≈ mpm_leslie.F
        @test lambda(mpm_sparse) ≈ lambda(mpm_leslie)
    end

end
