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

include("Results/DataFrames.jl")

end # module Results
