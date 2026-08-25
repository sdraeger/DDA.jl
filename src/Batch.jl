"""Batch processing and result aggregation for DDA."""
module Batch

using Statistics
using ..Results: STResult, CTResult, DEResult
import ..Results: n_channels, n_windows, n_coeffs, n_pairs
using ..API: run_st, run_ct, run_de
using ..DDADefaults
using ..Runner

# =============================================================================
# GroupResult
# =============================================================================

"""
    GroupResult

Aggregated DDA coefficient results from multiple subjects/recordings.

# Fields
- `coefficients::Array{Float64,4}`: Shape `(n_subjects, n_channels, n_windows, n_coeffs)`.
- `errors::Array{Float64,3}`: Shape `(n_subjects, n_channels, n_windows)`.
- `subject_labels::Vector{String}`: Label per subject.
- `channel_labels::Vector{String}`: Label per channel (or pair).
- `params::Dict{String,Any}`: Analysis parameters.
"""
struct GroupResult
    coefficients::Array{Float64,4}
    errors::Array{Float64,3}
    subject_labels::Vector{String}
    channel_labels::Vector{String}
    params::Dict{String,Any}
end

"""
    DEGroupResult

Aggregated dynamical-ergodicity results from multiple subjects/recordings,
returned by [`collect_results`](@ref) when collecting `DEResult`s.

# Fields
- `ergodicity::Matrix{Float64}`: Shape `(n_subjects, n_windows)`.
- `subject_labels::Vector{String}`: Label per subject.
- `params::Dict{String,Any}`: Analysis parameters.
"""
struct DEGroupResult
    ergodicity::Matrix{Float64}
    subject_labels::Vector{String}
    params::Dict{String,Any}
end

"""Number of subjects."""
n_subjects(g::GroupResult)::Int = size(g.coefficients, 1)

"""Number of subjects."""
n_subjects(g::DEGroupResult)::Int = size(g.ergodicity, 1)

"""Number of channels (or pairs)."""
n_channels(g::GroupResult)::Int = size(g.coefficients, 2)

"""Number of windows (truncated to minimum across subjects)."""
n_windows(g::GroupResult)::Int = size(g.coefficients, 3)

"""Number of windows (truncated to minimum across subjects)."""
n_windows(g::DEGroupResult)::Int = size(g.ergodicity, 2)

"""Number of coefficients."""
n_coeffs(g::GroupResult)::Int = size(g.coefficients, 4)

"""
    mean_over_windows(g::GroupResult) -> Array{Float64,3}

Average coefficients across the window dimension.

Returns shape `(n_subjects, n_channels, n_coeffs)`.
"""
function mean_over_windows(g::GroupResult)::Array{Float64,3}
    return dropdims(mean(g.coefficients; dims=3); dims=3)
end

"""
    mean_over_windows(g::DEGroupResult) -> Vector{Float64}

Average ergodicity per subject across the window dimension.

Returns shape `(n_subjects,)`.
"""
function mean_over_windows(g::DEGroupResult)::Vector{Float64}
    return vec(mean(g.ergodicity; dims=2))
end

# =============================================================================
# run_batch
# =============================================================================

"""
    run_batch(files; variant="st", kwargs...) -> Vector{Union{STResult, CTResult, DEResult}}

Process multiple ASCII files through DDA.

# Arguments
- `files::AbstractVector{<:AbstractString}`: Paths to ASCII data files.
- `variant::String="st"`: Analysis variant ("st", "ct", or "de").
- `delays`, `model`, `WL`, `WS`: Standard DDA parameters.
- `channel_labels`: Optional channel labels.
- `binary_path`: Explicit binary path.
- `progress::Bool=true`: Print progress.
- `load_func`: Custom file loader `f(path) -> Matrix{Float64}` (channels × samples).

# Returns
A list of result objects, one per file.
"""
function run_batch(
    files::AbstractVector{<:AbstractString};
    variant::String="st",
    delays::Vector{Int}=collect(DDADefaults.DELAYS),
    model::Runner.OptionalModelSpec=nothing,
    WL::Union{Int,Nothing}=DDADefaults.WL,
    WS::Union{Int,Nothing}=DDADefaults.WS,
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
        :WL => WL,
        :WS => WS,
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

    results = Vector{Any}(undef, length(files))
    parallel = Threads.nthreads() > 1 && length(files) > 1

    process = function(i)
        filepath = files[i]
        progress && println("[$i/$(length(files))] Processing $filepath")
        data = if load_func !== nothing
            load_func(filepath)          # assumed thread-safe when provided
        else
            _load_ascii(filepath)
        end
        results[i] = run_fn(; data=data, run_kwargs...)
    end

    if parallel
        Threads.@threads :dynamic for i in eachindex(files)
            process(i)
        end
    else
        foreach(process, eachindex(files))
    end

    return results
end

"""Load an ASCII file as Matrix{Float64} (channels × samples)."""
function _load_ascii(filepath::String)::Matrix{Float64}
    # Read tab/space-delimited numeric data, transpose to channels × samples
    rows = Runner._read_numeric_rows(filepath)
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

# =============================================================================
# collect_results
# =============================================================================

"""
    collect_results(results; labels=nothing) -> Union{GroupResult, DEGroupResult}

Stack multiple DDA results into a single aggregate. `STResult`s and `CTResult`s
produce a [`GroupResult`](@ref) with a 4D coefficient array; `DEResult`s produce
a [`DEGroupResult`](@ref).

Windows are truncated to the minimum count across subjects.

# Arguments
- `results`: List of `STResult`, `CTResult`, or `DEResult` (all same type).
- `labels::Union{Vector{String},Nothing}`: Subject labels (defaults to "subj_1", ...).
"""
function collect_results(
    results::Vector;
    labels::Union{Vector{String},Nothing}=nothing,
)::Union{GroupResult, DEGroupResult}
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
    return _collect_coeff_results(
        STResult,
        results,
        labels,
        n_channels,
        r -> r.channel_labels,
        "Channel",
    )
end

function _collect_ct(results::Vector, labels::Vector{String})::GroupResult
    return _collect_coeff_results(
        CTResult,
        results,
        labels,
        n_pairs,
        r -> r.pair_labels,
        "Pair",
    )
end

function _collect_coeff_results(
    ::Type{T},
    results::Vector,
    labels::Vector{String},
    entity_count::Function,
    entity_labels::Function,
    entity_name::String,
)::GroupResult where {T}
    n_subj = length(results)
    r1 = results[1]::T
    n_entities = entity_count(r1)
    nk = n_coeffs(r1)
    min_win = minimum(n_windows(r)::Int for r in results)

    for r in results
        r_typed = r::T
        entity_count(r_typed) == n_entities || error("$entity_name count mismatch")
        n_coeffs(r_typed) == nk || error("Coefficient count mismatch")
    end

    if any(n_windows(r)::Int != min_win for r in results)
        @warn "Window counts vary across subjects; truncating to $min_win"
    end

    coeffs = Array{Float64,4}(undef, n_subj, n_entities, min_win, nk)
    errs = Array{Float64,3}(undef, n_subj, n_entities, min_win)

    for (si, r) in enumerate(results)
        r_typed = r::T
        coeffs[si, :, :, :] = r_typed.coefficients[:, 1:min_win, :]
        errs[si, :, :] = r_typed.errors[:, 1:min_win]
    end

    return GroupResult(coeffs, errs, labels, entity_labels(r1), r1.params)
end

function _collect_de(results::Vector, labels::Vector{String})::DEGroupResult
    n_subj = length(results)
    min_win = minimum(n_windows(r)::Int for r in results)

    if any(n_windows(r)::Int != min_win for r in results)
        @warn "Window counts vary across subjects; truncating to $min_win"
    end

    ergodicity = Matrix{Float64}(undef, n_subj, min_win)
    for (si, r) in enumerate(results)
        ergodicity[si, :] = (r::DEResult).ergodicity[1:min_win]
    end

    return DEGroupResult(ergodicity, labels, (results[1]::DEResult).params)
end

end # module Batch
