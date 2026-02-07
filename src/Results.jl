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
- `window_starts::Vector{Int64}`: Start sample for each window.
- `window_ends::Vector{Int64}`: End sample for each window.
- `channel_labels::Vector{String}`: Label per channel.
- `params::Dict{String,Any}`: Analysis parameters used.
"""
struct STResult
    coefficients::Array{Float64,3}
    errors::Matrix{Float64}
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
- `window_starts::Vector{Int64}`: Start sample for each window.
- `window_ends::Vector{Int64}`: End sample for each window.
- `pair_labels::Vector{String}`: Label per channel pair.
- `params::Dict{String,Any}`: Analysis parameters used.
"""
struct CTResult
    coefficients::Array{Float64,3}
    errors::Matrix{Float64}
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
- `window_starts::Vector{Int64}`: Start sample for each window.
- `window_ends::Vector{Int64}`: End sample for each window.
- `params::Dict{String,Any}`: Analysis parameters used.
"""
struct DEResult
    ergodicity::Vector{Float64}
    window_starts::Vector{Int64}
    window_ends::Vector{Int64}
    params::Dict{String,Any}
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
    DF = @eval DataFrames
    nch = n_channels(r)
    nw = n_windows(r)
    nc = n_coeffs(r)

    rows = nch * nw
    col_channel = Vector{String}(undef, rows)
    col_wstart = Vector{Int64}(undef, rows)
    col_wend = Vector{Int64}(undef, rows)
    coeff_cols = [Vector{Float64}(undef, rows) for _ in 1:nc]
    col_error = Vector{Float64}(undef, rows)

    idx = 0
    for ch in 1:nch, w in 1:nw
        idx += 1
        col_channel[idx] = r.channel_labels[ch]
        col_wstart[idx] = r.window_starts[w]
        col_wend[idx] = r.window_ends[w]
        for c in 1:nc
            coeff_cols[c][idx] = r.coefficients[ch, w, c]
        end
        col_error[idx] = r.errors[ch, w]
    end

    pairs = [Symbol("channel") => col_channel,
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
    DF = @eval DataFrames
    np = n_pairs(r)
    nw = n_windows(r)
    nc = n_coeffs(r)

    rows = np * nw
    col_pair = Vector{String}(undef, rows)
    col_wstart = Vector{Int64}(undef, rows)
    col_wend = Vector{Int64}(undef, rows)
    coeff_cols = [Vector{Float64}(undef, rows) for _ in 1:nc]
    col_error = Vector{Float64}(undef, rows)

    idx = 0
    for p in 1:np, w in 1:nw
        idx += 1
        col_pair[idx] = r.pair_labels[p]
        col_wstart[idx] = r.window_starts[w]
        col_wend[idx] = r.window_ends[w]
        for c in 1:nc
            coeff_cols[c][idx] = r.coefficients[p, w, c]
        end
        col_error[idx] = r.errors[p, w]
    end

    pairs = [Symbol("pair") => col_pair,
             Symbol("window_start") => col_wstart,
             Symbol("window_end") => col_wend]
    for c in 1:nc
        push!(pairs, Symbol("a_$c") => coeff_cols[c])
    end
    push!(pairs, Symbol("error") => col_error)

    return Base.invokelatest(DF.DataFrame, pairs...)
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

    col_wstart = r.window_starts
    col_wend = r.window_ends
    col_erg = r.ergodicity

    return Base.invokelatest(
        DF.DataFrame,
        Symbol("window_start") => col_wstart,
        Symbol("window_end") => col_wend,
        Symbol("ergodicity") => col_erg,
    )
end

end # module Results
