"""Statistical testing for DDA results: permutation tests, effect sizes, window comparisons."""
module Stats

using Statistics
using Random
using HypothesisTests: UnequalVarianceTTest, MannWhitneyUTest, pvalue
using ..Results: STResult, CTResult, n_channels, n_pairs, n_coeffs

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
# Built-in statistical tests
# =============================================================================

"""Welch's unequal-variance t-test. Returns `(t_stat, p_value)`."""
function _welch_ttest(x::Vector{Float64}, y::Vector{Float64})::Tuple{Float64,Float64}
    nx, ny = length(x), length(y)
    (nx < 2 || ny < 2) && return (0.0, 1.0)

    mx, my = mean(x), mean(y)
    vx = var(x; corrected=true)
    vy = var(y; corrected=true)

    se = sqrt(vx / nx + vy / ny)
    se == 0.0 && return (0.0, 1.0)

    t_stat = (mx - my) / se
    return t_stat, pvalue(UnequalVarianceTTest(x, y))
end

"""Wilcoxon rank-sum (Mann-Whitney U). Returns `(U_stat, p_value)`."""
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

    u_stat = sum(ranks[1:nx]) - nx * (nx + 1) / 2

    p = try
        pvalue(MannWhitneyUTest(x, y))
    catch
        1.0
    end
    return (u_stat, p)
end

end # module Stats
