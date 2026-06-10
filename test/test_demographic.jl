@testset "Demographic stochasticity" begin
    rng = Random.Xoshiro(424242)

    # Two-stage model with explicit U (survival/growth) and F (fecundity).
    U = [0.0  0.0;
         0.6  0.55]            # juvenile -> adult (0.6); adult persists (0.55)
    F = [0.0  2.0;             # adults produce juveniles
         0.0  0.0]
    mpm = MatrixProjectionModel(U, F)
    A = mpm.A                  # == U + F
    n0 = [50, 50]

    @testset "ensemble mean tracks A^t n0 (multitype Galton-Watson)" begin
        prob = MPMProblem(Demographic(), mpm, n0, (0, 6))
        reps = 4000
        acc = [zeros(2) for _ in 1:7]
        for _ in 1:reps
            sol = solve(prob, DirectIteration(); rng=rng)
            for tt in 1:7
                acc[tt] .+= sol.u[tt]
            end
        end
        emp = [a ./ reps for a in acc]

        det = Vector{Vector{Float64}}(undef, 7)
        det[1] = Float64.(n0)
        for tt in 2:7
            det[tt] = A * det[tt-1]
        end

        for tt in 1:7
            @test isapprox(emp[tt], det[tt]; rtol=0.07)
        end
    end

    @testset "integer counts and genuine between-realization variance" begin
        prob = MPMProblem(Demographic(), mpm, n0, (0, 5))
        s1 = solve(prob, DirectIteration(); rng=Random.Xoshiro(1))
        s2 = solve(prob, DirectIteration(); rng=Random.Xoshiro(2))
        @test all(x -> x == round(x) && x >= 0, s1.u[end])   # integer-valued counts
        @test s1.u[end] != s2.u[end]                          # realizations differ
        @test length(s1.t) == 6 && length(s1.lambdas) == 5
    end

    @testset "subcritical model goes extinct (ensemble + quasi_extinction)" begin
        Usub = [0.0 0.0; 0.3 0.2]
        Fsub = [0.0 0.4; 0.0 0.0]
        mpm_sub = MatrixProjectionModel(Usub, Fsub)
        @test lambda(mpm_sub.A) < 1                           # subcritical
        prob = MPMProblem(Demographic(), mpm_sub, [20, 20], (0, 80))
        totals, sols = demographic_ensemble(prob; n_reps=400, rng=rng)
        @test size(totals) == (81, 400)
        qe = quasi_extinction(totals; threshold=1.0)
        @test qe.prob_extinct > 0.5
    end

    @testset "requires an A = U + F + C decomposition" begin
        # bare matrix (no decomposition) is rejected
        bareprob = MPMProblem(Demographic(), [0.0 2.0; 0.6 0.55], n0, (0, 3))
        @test_throws ErrorException solve(bareprob, DirectIteration())
        # fecundity folded into survival (column sum > 1) is rejected
        badmpm = MatrixProjectionModel([0.0 2.0; 0.6 0.55])   # U = A, col 2 sums to 2.55
        badprob = MPMProblem(Demographic(), badmpm, n0, (0, 3))
        @test_throws ErrorException solve(badprob, DirectIteration())
    end
end
