"""Statistical testing for DDA results: permutation tests, effect sizes, window comparisons."""
module Stats

using Statistics
using Random
using ..Results

export PermutationResult, EffectSizeResult, WindowComparisonResult
export permutation_test, compute_effect_size, compare_windows

# =============================================================================
# Result Types
# =============================================================================

"""
    PermutationResult

Result of a permutation test between two groups.

# Fields
- `observed_stat::Matrix{Float64}`: Shape `(n_channels, n_coeffs)`.
- `p_value::Matrix{Float64}`: Shape `(n_channels, n_coeffs)`.
- `null_distribution::Array{Float64,3}`: Shape `(n_permutations, n_channels, n_coeffs)`.
- `n_permutations::Int`: Number of permutations.
- `tail::Int`: 0=two-sided, 1=greater, -1=less.
"""
struct PermutationResult
    observed_stat::Matrix{Float64}
    p_value::Matrix{Float64}
    null_distribution::Array{Float64,3}
    n_permutations::Int
    tail::Int
end

"""
    EffectSizeResult

Cohen's d effect size between two groups.

# Fields
- `cohens_d::Matrix{Float64}`: Shape `(n_channels, n_coeffs)`.
- `mean_a::Matrix{Float64}`: Group A mean.
- `mean_b::Matrix{Float64}`: Group B mean.
"""
struct EffectSizeResult
    cohens_d::Matrix{Float64}
    mean_a::Matrix{Float64}
    mean_b::Matrix{Float64}
end

"""
    WindowComparisonResult

Within-subject comparison of baseline vs test windows.

# Fields
- `statistic::Matrix{Float64}`: Test statistic per (channel, coeff).
- `p_value::Matrix{Float64}`: P-value per (channel, coeff).
- `baseline_mean::Matrix{Float64}`: Mean of baseline windows.
- `test_mean::Matrix{Float64}`: Mean of test windows.
- `method::String`: "ttest" or "ranksum".
"""
struct WindowComparisonResult
    statistic::Matrix{Float64}
    p_value::Matrix{Float64}
    baseline_mean::Matrix{Float64}
    test_mean::Matrix{Float64}
    method::String
end

# =============================================================================
# Helpers
# =============================================================================

"""Extract window-averaged coefficients from a list of results.

Returns `(data, n_channels, n_coeffs)` where data is shape `(n_subjects, n_ch, n_coeff)`.
"""
_n_entities(r::STResult) = n_channels(r)
_n_entities(r::CTResult) = n_pairs(r)

function _extract_coefficients(results::Vector)::Array{Float64,3}
    isempty(results) && error("results must be non-empty")
    r1 = results[1]
    (r1 isa STResult || r1 isa CTResult) || error("permutation_test/compute_effect_size require STResult or CTResult")

    n_entities = _n_entities(r1)
    nk = n_coeffs(r1)
    data = Array{Float64,3}(undef, length(results), n_entities, nk)

    for (si, r) in enumerate(results)
        r_typed = r::typeof(r1)
        _n_entities(r_typed) == n_entities || error("Result entity count mismatch")
        n_coeffs(r_typed) == nk || error("Coefficient count mismatch")
        data[si, :, :] = dropdims(mean(r_typed.coefficients; dims=2); dims=2)
    end

    return data
end

"""Default test statistic: difference of means."""
function _default_stat(a::Array{Float64,3}, b::Array{Float64,3})::Matrix{Float64}
    return dropdims(mean(a; dims=1); dims=1) .- dropdims(mean(b; dims=1); dims=1)
end

# =============================================================================
# permutation_test
# =============================================================================

"""
    permutation_test(group_a, group_b; n_permutations=10000, stat_fun=nothing,
                     tail=0, seed=nothing) -> PermutationResult

Non-parametric permutation test comparing two groups of DDA results.

# Arguments
- `group_a`, `group_b`: Lists of `STResult` or `CTResult`.
- `n_permutations::Int=10000`: Number of permutations.
- `stat_fun`: Custom statistic function `(a, b) -> Matrix{Float64}`.
- `tail::Int=0`: 0=two-sided, 1=greater, -1=less.
- `seed`: Random seed for reproducibility.

# Returns
A [`PermutationResult`](@ref) with observed statistic, p-values, and null distribution.
"""
function permutation_test(
    group_a::Vector,
    group_b::Vector;
    n_permutations::Int=10000,
    stat_fun::Union{Function,Nothing}=nothing,
    tail::Int=0,
    seed::Union{Int,Nothing}=nothing,
)::PermutationResult
    tail in (-1, 0, 1) || error("tail must be -1, 0, or 1")

    data_a = _extract_coefficients(group_a)
    data_b = _extract_coefficients(group_b)

    _, nc, nk = size(data_a)
    size(data_a)[2:3] == size(data_b)[2:3] || error("Groups have different channel/coeff dimensions")

    sfun = something(stat_fun, _default_stat)

    observed = sfun(data_a, data_b)
    n_a = size(data_a, 1)
    n_total = n_a + size(data_b, 1)

    combined = cat(data_a, data_b; dims=1)
    null_dist = Array{Float64,3}(undef, n_permutations, nc, nk)

    rng = seed === nothing ? Random.default_rng() : MersenneTwister(seed)

    for p in 1:n_permutations
        perm = randperm(rng, n_total)
        perm_a = combined[perm[1:n_a], :, :]
        perm_b = combined[perm[(n_a+1):end], :, :]
        null_dist[p, :, :] = sfun(perm_a, perm_b)
    end

    p_value = Matrix{Float64}(undef, nc, nk)
    for ci in 1:nc, ki in 1:nk
        obs = observed[ci, ki]
        null_vals = null_dist[:, ci, ki]
        if tail == 0
            count = sum(abs.(null_vals) .>= abs(obs))
        elseif tail == 1
            count = sum(null_vals .>= obs)
        else  # tail == -1
            count = sum(null_vals .<= obs)
        end
        p_value[ci, ki] = (count + 1) / (n_permutations + 1)
    end

    return PermutationResult(observed, p_value, null_dist, n_permutations, tail)
end

# =============================================================================
# compute_effect_size
# =============================================================================

"""
    compute_effect_size(group_a, group_b) -> EffectSizeResult

Compute Cohen's d effect size between two groups of DDA results.

# Arguments
- `group_a`, `group_b`: Lists of `STResult` or `CTResult`.

# Returns
An [`EffectSizeResult`](@ref) with Cohen's d per (channel, coefficient).
"""
function compute_effect_size(
    group_a::Vector,
    group_b::Vector,
)::EffectSizeResult
    data_a = _extract_coefficients(group_a)
    data_b = _extract_coefficients(group_b)

    size(data_a)[2:3] == size(data_b)[2:3] || error("Groups have different dimensions")

    mean_a = dropdims(mean(data_a; dims=1); dims=1)
    mean_b = dropdims(mean(data_b; dims=1); dims=1)

    n_a = size(data_a, 1)
    n_b = size(data_b, 1)

    var_a = dropdims(var(data_a; dims=1, corrected=true); dims=1)
    var_b = dropdims(var(data_b; dims=1, corrected=true); dims=1)

    pooled_std = sqrt.(((n_a - 1) .* var_a .+ (n_b - 1) .* var_b) ./ (n_a + n_b - 2))

    cohens_d = similar(mean_a)
    for i in eachindex(pooled_std)
        if pooled_std[i] == 0.0
            cohens_d[i] = 0.0
        else
            cohens_d[i] = (mean_a[i] - mean_b[i]) / pooled_std[i]
        end
    end

    return EffectSizeResult(cohens_d, mean_a, mean_b)
end

# =============================================================================
# compare_windows
# =============================================================================

"""
    compare_windows(result, baseline_windows, test_windows; method="ttest") -> WindowComparisonResult

Compare baseline vs test windows within a single recording.

# Arguments
- `result::Union{STResult, CTResult}`: Single recording result.
- `baseline_windows`: Indices or `UnitRange` of baseline windows.
- `test_windows`: Indices or `UnitRange` of test windows.
- `method::String="ttest"`: `"ttest"` (Welch's t-test) or `"ranksum"` (Wilcoxon rank-sum).

# Returns
A [`WindowComparisonResult`](@ref) with test statistics and p-values.
"""
function compare_windows(
    result::Union{STResult,CTResult},
    baseline_windows,
    test_windows;
    method::String="ttest",
)::WindowComparisonResult
    method in ("ttest", "ranksum") || error("method must be 'ttest' or 'ranksum', got '$method'")

    coeffs = result.coefficients
    nc = size(coeffs, 1)
    nk = size(coeffs, 3)

    baseline = coeffs[:, baseline_windows, :]
    test_data = coeffs[:, test_windows, :]

    stat_mat = Matrix{Float64}(undef, nc, nk)
    p_mat = Matrix{Float64}(undef, nc, nk)

    for ci in 1:nc, ki in 1:nk
        b_vals = vec(baseline[ci, :, ki])
        t_vals = vec(test_data[ci, :, ki])

        if method == "ttest"
            t_stat, p_val = _welch_ttest(b_vals, t_vals)
        else
            t_stat, p_val = _ranksum(b_vals, t_vals)
        end
        stat_mat[ci, ki] = t_stat
        p_mat[ci, ki] = p_val
    end

    b_mean = dropdims(mean(baseline; dims=2); dims=2)
    t_mean = dropdims(mean(test_data; dims=2); dims=2)

    return WindowComparisonResult(stat_mat, p_mat, b_mean, t_mean, method)
end

# =============================================================================
# Built-in statistical tests (no external deps)
# =============================================================================

"""Welch's t-test for unequal variances. Returns (t_stat, p_value)."""
function _welch_ttest(x::Vector{Float64}, y::Vector{Float64})::Tuple{Float64,Float64}
    nx, ny = length(x), length(y)
    (nx < 2 || ny < 2) && return (0.0, 1.0)

    mx, my = mean(x), mean(y)
    vx = var(x; corrected=true)
    vy = var(y; corrected=true)

    se = sqrt(vx / nx + vy / ny)
    se == 0.0 && return (0.0, 1.0)

    t_stat = (mx - my) / se

    # Welch–Satterthwaite degrees of freedom
    num = (vx / nx + vy / ny)^2
    denom = (vx / nx)^2 / (nx - 1) + (vy / ny)^2 / (ny - 1)
    denom == 0.0 && return (t_stat, 1.0)
    df = num / denom

    # Approximate two-sided p-value using Student's t CDF approximation
    p = _t_pvalue(abs(t_stat), df)
    return (t_stat, p)
end

"""Wilcoxon rank-sum test (Mann-Whitney U). Returns (U_stat, p_value)."""
function _ranksum(x::Vector{Float64}, y::Vector{Float64})::Tuple{Float64,Float64}
    nx, ny = length(x), length(y)
    (nx == 0 || ny == 0) && return (0.0, 1.0)

    combined = vcat(x, y)
    n = nx + ny
    order = sortperm(combined)
    ranks = Vector{Float64}(undef, n)

    i = 1
    while i <= n
        j = i
        while j < n && combined[order[j+1]] == combined[order[j]]
            j += 1
        end
        avg_rank = (i + j) / 2.0
        for k in i:j
            ranks[order[k]] = avg_rank
        end
        i = j + 1
    end

    rank_sum_x = sum(ranks[1:nx])
    u_x = rank_sum_x - nx * (nx + 1) / 2

    # Normal approximation for large samples
    mu = nx * ny / 2.0
    sigma = sqrt(nx * ny * (n + 1) / 12.0)
    sigma == 0.0 && return (u_x, 1.0)

    z = (u_x - mu) / sigma
    p = 2.0 * _normal_cdf(-abs(z))
    return (u_x, p)
end

"""Approximate two-sided p-value for Student's t distribution."""
function _t_pvalue(t::Float64, df::Float64)::Float64
    # Use normal approximation for large df, otherwise Abramowitz & Stegun approx
    if df > 100
        return 2.0 * _normal_cdf(-t)
    end
    x = df / (df + t^2)
    # Regularized incomplete beta function approximation
    p = _betainc(df / 2, 0.5, x)
    return p
end

"""Standard normal CDF using Abramowitz & Stegun approximation (max error ~1.5e-7)."""
function _normal_cdf(x::Float64)::Float64
    # Rational approximation for Phi(x)
    if x < -8.0
        return 0.0
    elseif x > 8.0
        return 1.0
    end
    # Use symmetry for negative x
    if x < 0.0
        return 1.0 - _normal_cdf(-x)
    end
    # Abramowitz & Stegun 26.2.17
    b0 = 0.2316419
    b1 = 0.319381530
    b2 = -0.356563782
    b3 = 1.781477937
    b4 = -1.821255978
    b5 = 1.330274429
    t = 1.0 / (1.0 + b0 * x)
    phi = exp(-0.5 * x * x) / sqrt(2.0 * pi)
    return 1.0 - phi * t * (b1 + t * (b2 + t * (b3 + t * (b4 + t * b5))))
end

"""Log-gamma via Stirling's approximation (Lanczos, g=7, n=9). No external deps."""
function _lgamma(x::Float64)::Float64
    if x <= 0.0
        return Inf
    end
    # Lanczos coefficients (g=7, n=9)
    c = (
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7,
    )
    if x < 0.5
        return log(pi / sin(pi * x)) - _lgamma(1.0 - x)
    end
    x -= 1.0
    ag = c[1]
    for i in 2:9
        ag += c[i] / (x + i - 1)
    end
    t = x + 7.5
    return 0.5 * log(2.0 * pi) + (x + 0.5) * log(t) - t + log(ag)
end

"""
Regularized incomplete beta function I_x(a, b) via continued fraction.
Used for computing t-distribution p-values.
"""
function _betainc(a::Float64, b::Float64, x::Float64)::Float64
    (x <= 0.0) && return 0.0
    (x >= 1.0) && return 1.0

    # Use symmetry relation when x > (a+1)/(a+b+2)
    if x > (a + 1) / (a + b + 2)
        return 1.0 - _betainc(b, a, 1.0 - x)
    end

    # Lentz's continued fraction
    lbeta = _lgamma(a) + _lgamma(b) - _lgamma(a + b)
    front = exp(a * log(x) + b * log(1.0 - x) - lbeta) / a

    f = 1.0
    c = 1.0
    d = 1.0 - (a + b) * x / (a + 1)
    if abs(d) < 1e-30
        d = 1e-30
    end
    d = 1.0 / d
    f = d

    for m in 1:200
        # Even step
        num = m * (b - m) * x / ((a + 2m - 1) * (a + 2m))
        d = 1.0 + num * d
        abs(d) < 1e-30 && (d = 1e-30)
        c = 1.0 + num / c
        abs(c) < 1e-30 && (c = 1e-30)
        d = 1.0 / d
        f *= c * d

        # Odd step
        num = -(a + m) * (a + b + m) * x / ((a + 2m) * (a + 2m + 1))
        d = 1.0 + num * d
        abs(d) < 1e-30 && (d = 1e-30)
        c = 1.0 + num / c
        abs(c) < 1e-30 && (c = 1e-30)
        d = 1.0 / d
        delta = c * d
        f *= delta

        abs(delta - 1.0) < 1e-10 && break
    end

    return front * f
end

end # module Stats
