@testset "Life Tables" begin
    # 3-stage MPM (juvenile, adult, senescent)
    U = [0.0 0.0 0.0; 0.5 0.3 0.0; 0.0 0.4 0.2]
    F = [0.0 0.5 3.0; 0.0 0.0 0.0; 0.0 0.0 0.0]

    @testset "mpm_to_lx" begin
        lx = mpm_to_lx(U)
        @test lx[1] ≈ 1.0
        @test all(diff(lx) .<= 0)  # Non-increasing
        @test lx[end] < 0.01 + 0.01  # Eventually truncated
    end

    @testset "mpm_to_px" begin
        px = mpm_to_px(U)
        @test all(0 .<= px .<= 1.0 .+ 1e-10)
        @test px[end] == 0.0  # Absorbing at end
    end

    @testset "mpm_to_hx" begin
        hx = mpm_to_hx(U)
        @test all(hx .>= 0)
    end

    @testset "mpm_to_mx" begin
        mx = mpm_to_mx(U, F)
        @test length(mx) == length(mpm_to_lx(U))
        @test mx[1] == 0.0  # No reproduction at age 0 (start in stage 1, no fecundity)
    end

    @testset "mpm_to_table" begin
        tbl = mpm_to_table(U, F)
        @test haskey(tbl, :x)
        @test haskey(tbl, :lx)
        @test haskey(tbl, :px)
        @test haskey(tbl, :hx)
        @test haskey(tbl, :mx)
        @test length(tbl.x) == length(tbl.lx)
    end

    @testset "Conversions roundtrip" begin
        lx = [1.0, 0.8, 0.5, 0.2, 0.05]

        # lx → px → lx
        px = lx_to_px(lx)
        lx2 = px_to_lx(px)
        @test lx2[1:length(lx)] ≈ lx atol=1e-10

        # lx → hx → lx
        hx = lx_to_hx(lx)
        lx3 = hx_to_lx(hx)
        @test lx3[1:length(lx)] ≈ lx atol=1e-10

        # px → hx → px
        hx2 = px_to_hx(px)
        px2 = hx_to_px(hx2)
        @test px2 ≈ px atol=1e-10
    end

    @testset "start parameter" begin
        lx_start1 = mpm_to_lx(U; start=1)
        lx_start2 = mpm_to_lx(U; start=2)
        @test lx_start1[1] == lx_start2[1]  # Both start at 1.0
        # Different starting stages should give different schedules
        @test lx_start1[3] != lx_start2[3]
    end
end
