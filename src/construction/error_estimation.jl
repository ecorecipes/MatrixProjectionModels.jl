"""
Error estimation via bootstrapping for MPM-derived statistics.
Based on mpmsim::calculate_errors and mpmsim::compute_ci.
"""

"""
    calculate_errors(mpm::MatrixProjectionModel, sample_size::Int;
                     type::Symbol=:sem, rng::AbstractRNG=Random.default_rng(),
                     n_boot::Int=1000)

Calculate standard errors or confidence intervals for MPM entries
via bootstrap resampling with sampling error.

# Arguments
- `type`: `:sem` for standard error of the mean, `:ci95` for 95% CI
- `n_boot`: Number of bootstrap replicates

# Returns
NamedTuple with (A, U, F, C) matrices of errors.
"""
function calculate_errors(mpm::MatrixProjectionModel, sample_size::Int;
                          type::Symbol=:sem,
                          rng::AbstractRNG=Random.default_rng(),
                          n_boot::Int=1000)
    n = n_stages(mpm)
    boot_A = zeros(n, n, n_boot)

    for b in 1:n_boot
        mpm_b = add_mpm_error(mpm, sample_size; rng=rng)
        boot_A[:, :, b] .= mpm_b.A
    end

    if type == :sem
        err_A = dropdims(std(boot_A; dims=3); dims=3)
        err_U = zeros(n, n)
        err_F = zeros(n, n)
        for i in 1:n, j in 1:n
            if mpm.U[i, j] > 0
                err_U[i, j] = err_A[i, j]
            end
            if mpm.F[i, j] > 0
                err_F[i, j] = err_A[i, j]
            end
        end
        return (A=err_A, U=err_U, F=err_F)
    elseif type == :ci95
        lo = dropdims(mapslices(x -> quantile(x, 0.025), boot_A; dims=3); dims=3)
        hi = dropdims(mapslices(x -> quantile(x, 0.975), boot_A; dims=3); dims=3)
        return (lower=lo, upper=hi)
    else
        throw(ArgumentError("type must be :sem or :ci95"))
    end
end

"""
    compute_ci(mpm::MatrixProjectionModel, fn::Function, sample_size::Int;
               n_boot::Int=1000, ci::Float64=0.95,
               rng::AbstractRNG=Random.default_rng())

Compute bootstrap confidence interval for any statistic `fn(mpm)`.

# Arguments
- `fn`: Function that takes a MatrixProjectionModel and returns a scalar
- `sample_size`: Number of individuals per stage for sampling error
- `n_boot`: Number of bootstrap replicates
- `ci`: Confidence level

# Returns
NamedTuple (estimate, lower, upper, se, boot_values).
"""
function compute_ci(mpm::MatrixProjectionModel, fn::Function, sample_size::Int;
                    n_boot::Int=1000, ci::Float64=0.95,
                    rng::AbstractRNG=Random.default_rng())
    estimate = fn(mpm)
    boot_vals = Float64[]

    for _ in 1:n_boot
        mpm_b = add_mpm_error(mpm, sample_size; rng=rng)
        try
            push!(boot_vals, fn(mpm_b))
        catch
            # Skip failed bootstrap samples
        end
    end

    alpha = (1 - ci) / 2
    lower = quantile(boot_vals, alpha)
    upper = quantile(boot_vals, 1 - alpha)
    se = std(boot_vals)

    return (estimate=estimate, lower=lower, upper=upper, se=se, boot_values=boot_vals)
end
