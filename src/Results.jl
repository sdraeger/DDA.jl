"""Structured result types for DDA analysis."""
module Results

export STResult, CTResult, DEResult
export n_channels, n_windows, n_coeffs, n_pairs, to_dataframe

"""
    STResult

Single Timeseries analysis result.

# Fields
- `coefficients::Array{Float64,3}`: Shape `(n_channels, n_windows, n_coeffs)`.
- `errors::Matrix{Float64}`: Shape `(n_channels, n_windows)`.
- `T::Vector{Float64}`: Raw first column from the DDA output file.
- `t::Vector{Float64}`: Derived time axis `(T + 1 + derivative_points + TM) / SR` when `sampling_rate` is provided, otherwise unscaled.
- `window_starts::Vector{Int64}`: Start sample for each window.
- `window_ends::Vector{Int64}`: End sample for each window.
- `channel_labels::Vector{String}`: Label per channel.
- `params::Dict{String,Any}`: Analysis parameters used.
"""
struct STResult
    coefficients::Array{Float64,3}
    errors::Matrix{Float64}
    T::Vector{Float64}
    t::Vector{Float64}
    window_starts::Vector{Int64}
    window_ends::Vector{Int64}
    channel_labels::Vector{String}
    params::Dict{String,Any}
end

"""
    CTResult

Cross-Timeseries analysis result.

# Fields
- `coefficients::Array{Float64,3}`: Shape `(n_pairs, n_windows, n_coeffs)`.
- `errors::Matrix{Float64}`: Shape `(n_pairs, n_windows)`.
- `T::Vector{Float64}`: Raw first column from the DDA output file.
- `t::Vector{Float64}`: Derived time axis `(T + 1 + derivative_points + TM) / SR` when `sampling_rate` is provided, otherwise unscaled.
- `window_starts::Vector{Int64}`: Start sample for each window.
- `window_ends::Vector{Int64}`: End sample for each window.
- `pair_labels::Vector{String}`: Label per channel pair.
- `params::Dict{String,Any}`: Analysis parameters used.
"""
struct CTResult
    coefficients::Array{Float64,3}
    errors::Matrix{Float64}
    T::Vector{Float64}
    t::Vector{Float64}
    window_starts::Vector{Int64}
    window_ends::Vector{Int64}
    pair_labels::Vector{String}
    params::Dict{String,Any}
end

"""
    DEResult

Dynamical Ergodicity analysis result.

# Fields
- `ergodicity::Vector{Float64}`: Ergodicity measure per window.
- `T::Vector{Float64}`: Raw first column from the DDA output file.
- `t::Vector{Float64}`: Derived time axis `(T + 1 + derivative_points + TM) / SR` when `sampling_rate` is provided, otherwise unscaled.
- `window_starts::Vector{Int64}`: Start sample for each window.
- `window_ends::Vector{Int64}`: End sample for each window.
- `params::Dict{String,Any}`: Analysis parameters used.
"""
struct DEResult
    ergodicity::Vector{Float64}
    T::Vector{Float64}
    t::Vector{Float64}
    window_starts::Vector{Int64}
    window_ends::Vector{Int64}
    params::Dict{String,Any}
end

function STResult(
    coefficients::Array{Float64,3},
    errors::Matrix{Float64},
    window_starts::Vector{<:Integer},
    window_ends::Vector{<:Integer},
    channel_labels::Vector{String},
    params::Dict{String,Any},
)
    raw_T = Float64.(window_starts)
    derived_t = Float64.(window_starts)
    return STResult(
        coefficients,
        errors,
        raw_T,
        derived_t,
        Int64.(window_starts),
        Int64.(window_ends),
        channel_labels,
        params,
    )
end

function CTResult(
    coefficients::Array{Float64,3},
    errors::Matrix{Float64},
    window_starts::Vector{<:Integer},
    window_ends::Vector{<:Integer},
    pair_labels::Vector{String},
    params::Dict{String,Any},
)
    raw_T = Float64.(window_starts)
    derived_t = Float64.(window_starts)
    return CTResult(
        coefficients,
        errors,
        raw_T,
        derived_t,
        Int64.(window_starts),
        Int64.(window_ends),
        pair_labels,
        params,
    )
end

function DEResult(
    ergodicity::Vector{Float64},
    window_starts::Vector{<:Integer},
    window_ends::Vector{<:Integer},
    params::Dict{String,Any},
)
    raw_T = Float64.(window_starts)
    derived_t = Float64.(window_starts)
    return DEResult(
        ergodicity,
        raw_T,
        derived_t,
        Int64.(window_starts),
        Int64.(window_ends),
        params,
    )
end

# ---------------------------------------------------------------------------
# Accessor functions (Julia multiple dispatch instead of Python @property)
# ---------------------------------------------------------------------------

"""Number of channels in an STResult."""
n_channels(r::STResult)::Int = size(r.coefficients, 1)

"""Number of time windows."""
n_windows(r::STResult)::Int = size(r.coefficients, 2)

"""Number of DDA coefficients per window."""
n_coeffs(r::STResult)::Int = size(r.coefficients, 3)

"""Number of channel pairs in a CTResult."""
n_pairs(r::CTResult)::Int = size(r.coefficients, 1)

"""Number of time windows."""
n_windows(r::CTResult)::Int = size(r.coefficients, 2)

"""Number of DDA coefficients per window."""
n_coeffs(r::CTResult)::Int = size(r.coefficients, 3)

"""Number of time windows."""
n_windows(r::DEResult)::Int = length(r.ergodicity)

# ---------------------------------------------------------------------------
# DataFrame conversion (lazy DataFrames.jl import)
# ---------------------------------------------------------------------------

# Module-level flag to track whether DataFrames has been loaded
const _dataframes_loaded = Ref(false)

function _ensure_dataframes()
    if !_dataframes_loaded[]
        try
            @eval import DataFrames
            _dataframes_loaded[] = true
        catch
            error(
                "DataFrames.jl is required for to_dataframe(). " *
                "Install with: using Pkg; Pkg.add(\"DataFrames\")"
            )
        end
    end
end

"""
    to_dataframe(r::STResult)

Convert an STResult to a long-format DataFrame.

Columns: `channel`, `window_start`, `window_end`, `a_1`, ..., `a_N`, `error`.

Requires DataFrames.jl to be installed.
"""
function to_dataframe(r::STResult)
    _ensure_dataframes()
    Base.invokelatest(_to_dataframe_st, r)
end

function _to_dataframe_st(r::STResult)
    return _to_dataframe_coeff_result(r, :channel, r.channel_labels)
end

function _to_dataframe_coeff_result(
    r::Union{STResult,CTResult},
    entity_column::Symbol,
    labels::Vector{String},
)
    DF = @eval DataFrames
    n_entities = size(r.coefficients, 1)
    nw = n_windows(r)
    nc = n_coeffs(r)

    rows = n_entities * nw
    col_entity = Vector{String}(undef, rows)
    col_T = Vector{Float64}(undef, rows)
    col_t = Vector{Float64}(undef, rows)
    col_wstart = Vector{Int64}(undef, rows)
    col_wend = Vector{Int64}(undef, rows)
    coeff_cols = [Vector{Float64}(undef, rows) for _ in 1:nc]
    col_error = Vector{Float64}(undef, rows)

    idx = 0
    for entity in 1:n_entities, w in 1:nw
        idx += 1
        col_entity[idx] = labels[entity]
        col_T[idx] = r.T[w]
        col_t[idx] = r.t[w]
        col_wstart[idx] = r.window_starts[w]
        col_wend[idx] = r.window_ends[w]
        for c in 1:nc
            coeff_cols[c][idx] = r.coefficients[entity, w, c]
        end
        col_error[idx] = r.errors[entity, w]
    end

    pairs = [entity_column => col_entity,
             Symbol("T") => col_T,
             Symbol("t") => col_t,
             Symbol("window_start") => col_wstart,
             Symbol("window_end") => col_wend]
    for c in 1:nc
        push!(pairs, Symbol("a_$c") => coeff_cols[c])
    end
    push!(pairs, Symbol("error") => col_error)

    return Base.invokelatest(DF.DataFrame, pairs...)
end

"""
    to_dataframe(r::CTResult)

Convert a CTResult to a long-format DataFrame.
"""
function to_dataframe(r::CTResult)
    _ensure_dataframes()
    Base.invokelatest(_to_dataframe_ct, r)
end

function _to_dataframe_ct(r::CTResult)
    return _to_dataframe_coeff_result(r, :pair, r.pair_labels)
end

"""
    to_dataframe(r::DEResult)

Convert a DEResult to a long-format DataFrame.
"""
function to_dataframe(r::DEResult)
    _ensure_dataframes()
    Base.invokelatest(_to_dataframe_de, r)
end

function _to_dataframe_de(r::DEResult)
    DF = @eval DataFrames
    nw = n_windows(r)

    col_T = r.T
    col_t = r.t
    col_wstart = r.window_starts
    col_wend = r.window_ends
    col_erg = r.ergodicity

    return Base.invokelatest(
        DF.DataFrame,
        Symbol("T") => col_T,
        Symbol("t") => col_t,
        Symbol("window_start") => col_wstart,
        Symbol("window_end") => col_wend,
        Symbol("ergodicity") => col_erg,
    )
end

end # module Results
