"""Batch processing and result aggregation for DDA."""
module Batch

using Statistics
using ..Results
using ..API

export GroupResult, run_batch, collect_results
export n_subjects, mean_over_windows

# =============================================================================
# GroupResult
# =============================================================================

"""
    GroupResult

Aggregated DDA results from multiple subjects/recordings.

# Fields
- `coefficients::Array{Float64,4}`: Shape `(n_subjects, n_channels, n_windows, n_coeffs)`.
- `errors::Array{Float64,3}`: Shape `(n_subjects, n_channels, n_windows)`.
- `subject_labels::Vector{String}`: Label per subject.
- `channel_labels::Vector{String}`: Label per channel (or pair).
- `params::Dict{String,Any}`: Analysis parameters.
- `variant::String`: Variant abbreviation ("ST", "CT", or "DE").
"""
struct GroupResult
    coefficients::Array{Float64,4}
    errors::Array{Float64,3}
    subject_labels::Vector{String}
    channel_labels::Vector{String}
    params::Dict{String,Any}
    variant::String
end

"""Number of subjects."""
n_subjects(g::GroupResult)::Int = size(g.coefficients, 1)

"""Number of channels (or pairs)."""
Results.n_channels(g::GroupResult)::Int = size(g.coefficients, 2)

"""Number of windows (truncated to minimum across subjects)."""
Results.n_windows(g::GroupResult)::Int = size(g.coefficients, 3)

"""Number of coefficients."""
Results.n_coeffs(g::GroupResult)::Int = size(g.coefficients, 4)

"""
    mean_over_windows(g::GroupResult) -> Array{Float64,3}

Average coefficients across the window dimension.

Returns shape `(n_subjects, n_channels, n_coeffs)`.
"""
function mean_over_windows(g::GroupResult)::Array{Float64,3}
    return dropdims(mean(g.coefficients; dims=3); dims=3)
end

# =============================================================================
# run_batch
# =============================================================================

"""
    run_batch(files; variant="st", kwargs...) -> Vector{Union{STResult, CTResult, DEResult}}

Process multiple ASCII files through DDA.

# Arguments
- `files::Vector{String}`: Paths to ASCII data files.
- `variant::String="st"`: Analysis variant ("st", "ct", or "de").
- `sfreq::Float64=1.0`: Sampling frequency.
- `delays`, `model`, `wl`, `ws`: Standard DDA parameters.
- `channel_labels`: Optional channel labels.
- `binary_path`: Explicit binary path.
- `progress::Bool=true`: Print progress.
- `load_func`: Custom file loader `f(path) -> Matrix{Float64}` (channels × samples).

# Returns
A list of result objects, one per file.
"""
function run_batch(
    files::Vector{String};
    variant::String="st",
    sfreq::Float64=1.0,
    delays::Vector{Int}=collect(Results.DDADefaults_DELAYS()),
    model::Union{Vector{Int},Nothing}=nothing,
    wl::Int=200,
    ws::Int=100,
    channel_labels::Union{Vector{String},Nothing}=nothing,
    binary_path::Union{String,Nothing}=nothing,
    progress::Bool=true,
    load_func::Union{Function,Nothing}=nothing,
    kwargs...,
)
    variant_lower = lowercase(variant)
    variant_lower in ("st", "ct", "de") || error("variant must be 'st', 'ct', or 'de', got '$variant'")

    for f in files
        isfile(f) || error("File not found: $f")
    end

    run_fn = variant_lower == "st" ? run_st :
             variant_lower == "ct" ? run_ct : run_de

    # Build kwargs for run function
    run_kwargs = Dict{Symbol,Any}(
        :sfreq => sfreq,
        :wl => wl,
        :ws => ws,
    )
    if model !== nothing
        run_kwargs[:model] = model
    end
    if !isempty(delays)
        run_kwargs[:delays] = delays
    end
    if channel_labels !== nothing
        run_kwargs[:channel_labels] = channel_labels
    end
    if binary_path !== nothing
        run_kwargs[:binary_path] = binary_path
    end
    # Pass through extra kwargs
    for (k, v) in kwargs
        run_kwargs[k] = v
    end

    results = []
    for (i, filepath) in enumerate(files)
        progress && println("[$i/$(length(files))] Processing $filepath")
        data = if load_func !== nothing
            load_func(filepath)
        else
            _load_ascii(filepath)
        end
        result = run_fn(data; run_kwargs...)
        push!(results, result)
    end

    return results
end

"""Load an ASCII file as Matrix{Float64} (channels × samples)."""
function _load_ascii(filepath::String)::Matrix{Float64}
    # Read tab/space-delimited numeric data, transpose to channels × samples
    lines = readlines(filepath)
    rows = Vector{Vector{Float64}}()
    for line in lines
        stripped = strip(line)
        (isempty(stripped) || startswith(stripped, '#')) && continue
        parts = split(stripped)
        values = tryparse.(Float64, parts)
        all(v -> v !== nothing, values) || continue
        push!(rows, Float64[v for v in values])
    end
    isempty(rows) && error("No numeric data in $filepath")
    n_timepoints = length(rows)
    n_channels = length(rows[1])
    mat = Matrix{Float64}(undef, n_channels, n_timepoints)
    for (t, row) in enumerate(rows)
        length(row) == n_channels || error("Inconsistent columns in $filepath at row $t")
        for ch in 1:n_channels
            mat[ch, t] = row[ch]
        end
    end
    return mat
end

# Workaround: access DDADefaults.DELAYS without a direct dependency
DDADefaults_DELAYS() = [7, 10]

# =============================================================================
# collect_results
# =============================================================================

"""
    collect_results(results; labels=nothing) -> GroupResult

Stack multiple DDA results into a single `GroupResult` with a 4D coefficient array.

Windows are truncated to the minimum count across subjects.

# Arguments
- `results`: List of `STResult`, `CTResult`, or `DEResult` (all same type).
- `labels::Union{Vector{String},Nothing}`: Subject labels (defaults to "subj_1", ...).

# Returns
A [`GroupResult`](@ref) with coefficients shape `(n_subjects, n_channels, min_windows, n_coeffs)`.
"""
function collect_results(
    results::Vector;
    labels::Union{Vector{String},Nothing}=nothing,
)::GroupResult
    isempty(results) && error("results must be non-empty")
    T = typeof(results[1])
    all(r -> typeof(r) == T, results) || error("All results must be the same type")

    n_subj = length(results)
    subj_labels = something(labels, ["subj_$i" for i in 1:n_subj])
    length(subj_labels) == n_subj || error("Number of labels must match number of results")

    if T === DEResult
        return _collect_de(results, subj_labels)
    elseif T === STResult
        return _collect_st(results, subj_labels)
    elseif T === CTResult
        return _collect_ct(results, subj_labels)
    else
        error("Unsupported result type: $T")
    end
end

function _collect_st(results::Vector, labels::Vector{String})::GroupResult
    n_subj = length(results)
    r1 = results[1]::STResult
    nc = n_channels(r1)
    nk = n_coeffs(r1)
    min_win = minimum(n_windows(r)::Int for r in results)

    for r in results
        r_typed = r::STResult
        n_channels(r_typed) == nc || error("Channel count mismatch: expected $nc, got $(n_channels(r_typed))")
        n_coeffs(r_typed) == nk || error("Coefficient count mismatch: expected $nk, got $(n_coeffs(r_typed))")
    end

    if any(n_windows(r)::Int != min_win for r in results)
        @warn "Window counts vary across subjects; truncating to $min_win"
    end

    coeffs = Array{Float64,4}(undef, n_subj, nc, min_win, nk)
    errs = Array{Float64,3}(undef, n_subj, nc, min_win)

    for (si, r) in enumerate(results)
        r_typed = r::STResult
        coeffs[si, :, :, :] = r_typed.coefficients[:, 1:min_win, :]
        errs[si, :, :] = r_typed.errors[:, 1:min_win]
    end

    ch_labels = r1.channel_labels
    params = r1.params
    return GroupResult(coeffs, errs, labels, ch_labels, params, "ST")
end

function _collect_ct(results::Vector, labels::Vector{String})::GroupResult
    n_subj = length(results)
    r1 = results[1]::CTResult
    np = n_pairs(r1)
    nk = n_coeffs(r1)
    min_win = minimum(n_windows(r)::Int for r in results)

    for r in results
        r_typed = r::CTResult
        n_pairs(r_typed) == np || error("Pair count mismatch")
        n_coeffs(r_typed) == nk || error("Coefficient count mismatch")
    end

    if any(n_windows(r)::Int != min_win for r in results)
        @warn "Window counts vary across subjects; truncating to $min_win"
    end

    coeffs = Array{Float64,4}(undef, n_subj, np, min_win, nk)
    errs = Array{Float64,3}(undef, n_subj, np, min_win)

    for (si, r) in enumerate(results)
        r_typed = r::CTResult
        coeffs[si, :, :, :] = r_typed.coefficients[:, 1:min_win, :]
        errs[si, :, :] = r_typed.errors[:, 1:min_win]
    end

    pair_labels = r1.pair_labels
    params = r1.params
    return GroupResult(coeffs, errs, labels, pair_labels, params, "CT")
end

function _collect_de(results::Vector, labels::Vector{String})::GroupResult
    n_subj = length(results)
    min_win = minimum(n_windows(r)::Int for r in results)

    if any(n_windows(r)::Int != min_win for r in results)
        @warn "Window counts vary across subjects; truncating to $min_win"
    end

    # DE: 1 "channel", 0 coefficients — store ergodicity in a dummy axis
    coeffs = Array{Float64,4}(undef, n_subj, 1, min_win, 1)
    errs = Array{Float64,3}(undef, n_subj, 1, min_win)

    for (si, r) in enumerate(results)
        r_typed = r::DEResult
        for w in 1:min_win
            coeffs[si, 1, w, 1] = r_typed.ergodicity[w]
            errs[si, 1, w] = 0.0
        end
    end

    r1 = results[1]::DEResult
    params = r1.params
    return GroupResult(coeffs, errs, labels, ["ergodicity"], params, "DE")
end

end # module Batch
