"""High-level analysis functions: run_st, run_ct, run_de."""
module API

using Printf
using ..DDADefaults
using ..Results
using ..Runner
using ..Variants

export run_st, run_ct, run_de

# =============================================================================
# HELPERS
# =============================================================================

"""Write a Matrix{Float64} (n_channels × n_samples) to a temporary ASCII file.

Returns the path to the written file. The file is transposed so each row is a
timepoint and each column is a channel (the format the DDA binary expects).
Uses fixed-point notation (no scientific 'e') since the DDA binary rejects letters.
"""
function _write_temp_ascii(data::Matrix{Float64})::String
    path = tempname() * ".txt"
    open(path, "w") do io
        n_ch, n_samp = size(data)
        for t in 1:n_samp
            for ch in 1:n_ch
                ch > 1 && print(io, '\t')
                @printf(io, "%.15f", data[ch, t])
            end
            println(io)
        end
    end
    return path
end

"""Build a params dict capturing analysis configuration."""
function _make_params(;
    sfreq::Float64,
    delays::Vector{Int},
    model::Vector{Int},
    wl::Int,
    ws::Int,
    dm::Int,
    order::Int,
    nr_tau::Int,
    extra::Dict{String,Any}=Dict{String,Any}(),
)::Dict{String,Any}
    d = Dict{String,Any}(
        "sfreq" => sfreq,
        "delays" => delays,
        "model" => model,
        "wl" => wl,
        "ws" => ws,
        "dm" => dm,
        "order" => order,
        "nr_tau" => nr_tau,
    )
    merge!(d, extra)
    return d
end

# =============================================================================
# run_st
# =============================================================================

"""
    run_st(data; sfreq=1.0, delays=..., kwargs...) -> STResult

Run Single Timeseries DDA on a `Matrix{Float64}` of shape `(n_channels, n_samples)`.

# Arguments
- `data::Matrix{Float64}`: Input data, channels × samples.
- `sfreq::Float64=1.0`: Sampling frequency in Hz.
- `delays::Vector{Int}`: Delay (tau) values.
- `model::Vector{Int}`: Model term indices.
- `wl::Int`: Window length.
- `ws::Int`: Window step.
- `channel_labels::Union{Vector{String},Nothing}`: Optional labels.
- `binary_path::Union{String,Nothing}`: Explicit binary path.
- `dm::Int`, `order::Int`, `nr_tau::Int`: Model parameters.

# Returns
An [`STResult`](@ref) with fields `coefficients`, `errors`, `window_starts`, etc.
"""
function run_st(
    data::Matrix{Float64};
    sfreq::Float64=1.0,
    delays::Vector{Int}=collect(DDADefaults.DELAYS),
    model::Vector{Int}=copy(DDADefaults.MODEL_PARAMS),
    wl::Int=DDADefaults.WINDOW_LENGTH,
    ws::Int=DDADefaults.WINDOW_STEP,
    channel_labels::Union{Vector{String},Nothing}=nothing,
    binary_path::Union{String,Nothing}=nothing,
    dm::Int=DDADefaults.MODEL_DIMENSION,
    order::Int=DDADefaults.POLYNOMIAL_ORDER,
    nr_tau::Int=DDADefaults.NUM_TAU,
)::STResult
    n_ch, n_samp = size(data)
    labels = something(channel_labels, ["ch$i" for i in 1:n_ch])

    tmp = _write_temp_ascii(data)
    try
        channels_0based = collect(0:(n_ch - 1))
        request = DDARequest(
            tmp, channels_0based, ["ST"];
            window_length=wl, window_step=ws,
            delays=delays, model_encoding=model,
            dm=dm, order=order, nr_tau=nr_tau,
            sampling_rate=sfreq > 1000.0 ? sfreq : nothing,
        )
        runner = binary_path === nothing ? DDARunner() : DDARunner(binary_path)
        raw = run_analysis_structured(runner, request)

        haskey(raw, "ST") || error("No ST results in DDA output")
        return _raw_to_st(raw["ST"], labels, sfreq, delays, model, wl, ws, dm, order, nr_tau)
    finally
        isfile(tmp) && rm(tmp; force=true)
    end
end

function _raw_to_st(
    channels::Vector{StructuredChannelData},
    labels::Vector{String},
    sfreq::Float64,
    delays::Vector{Int},
    model::Vector{Int},
    wl::Int, ws::Int, dm::Int, order::Int, nr_tau::Int,
)::STResult
    n_ch = length(channels)
    n_win = length(channels[1].timepoints)
    n_coeff = isempty(channels[1].timepoints) ? 0 : length(channels[1].timepoints[1].coefficients)

    coefficients = Array{Float64,3}(undef, n_ch, n_win, n_coeff)
    errors = Matrix{Float64}(undef, n_ch, n_win)
    win_starts = Vector{Int64}(undef, n_win)
    win_ends = Vector{Int64}(undef, n_win)

    for (ci, ch_data) in enumerate(channels)
        for (wi, tp) in enumerate(ch_data.timepoints)
            if ci == 1
                win_starts[wi] = tp.window_start
                win_ends[wi] = tp.window_end
            end
            for (ki, c) in enumerate(tp.coefficients)
                coefficients[ci, wi, ki] = c
            end
            errors[ci, wi] = tp.error
        end
    end

    params = _make_params(;
        sfreq=sfreq, delays=delays, model=model,
        wl=wl, ws=ws, dm=dm, order=order, nr_tau=nr_tau,
    )
    return STResult(coefficients, errors, win_starts, win_ends, labels, params)
end

# =============================================================================
# run_ct
# =============================================================================

"""
    run_ct(data; sfreq=1.0, delays=..., kwargs...) -> CTResult

Run Cross-Timeseries DDA. Requires at least 2 channels.

# Arguments
Same as `run_st`, plus:
- `ct_wl::Union{Int,Nothing}`: CT window length (defaults to `wl`).
- `ct_ws::Union{Int,Nothing}`: CT window step (defaults to `ws`).

# Returns
A [`CTResult`](@ref) with pair-wise analysis results.
"""
function run_ct(
    data::Matrix{Float64};
    sfreq::Float64=1.0,
    delays::Vector{Int}=collect(DDADefaults.DELAYS),
    model::Vector{Int}=copy(DDADefaults.MODEL_PARAMS),
    wl::Int=DDADefaults.WINDOW_LENGTH,
    ws::Int=DDADefaults.WINDOW_STEP,
    ct_wl::Union{Int,Nothing}=nothing,
    ct_ws::Union{Int,Nothing}=nothing,
    channel_labels::Union{Vector{String},Nothing}=nothing,
    binary_path::Union{String,Nothing}=nothing,
    dm::Int=DDADefaults.MODEL_DIMENSION,
    order::Int=DDADefaults.POLYNOMIAL_ORDER,
    nr_tau::Int=DDADefaults.NUM_TAU,
)::CTResult
    n_ch, _ = size(data)
    n_ch >= 2 || error("CT analysis requires at least 2 channels, got $n_ch")

    labels = something(channel_labels, ["ch$i" for i in 1:n_ch])
    pair_labels = String[]
    for i in 1:n_ch, j in (i+1):n_ch
        push!(pair_labels, "$(labels[i])-$(labels[j])")
    end

    tmp = _write_temp_ascii(data)
    try
        channels_0based = collect(0:(n_ch - 1))
        request = DDARequest(
            tmp, channels_0based, ["CT"];
            window_length=wl, window_step=ws,
            ct_window_length=something(ct_wl, wl),
            ct_window_step=something(ct_ws, ws),
            delays=delays, model_encoding=model,
            dm=dm, order=order, nr_tau=nr_tau,
            sampling_rate=sfreq > 1000.0 ? sfreq : nothing,
        )
        runner = binary_path === nothing ? DDARunner() : DDARunner(binary_path)
        raw = run_analysis_structured(runner, request)

        haskey(raw, "CT") || error("No CT results in DDA output")
        return _raw_to_ct(raw["CT"], pair_labels, sfreq, delays, model, wl, ws, dm, order, nr_tau)
    finally
        isfile(tmp) && rm(tmp; force=true)
    end
end

function _raw_to_ct(
    pairs::Vector{StructuredChannelData},
    pair_labels::Vector{String},
    sfreq::Float64,
    delays::Vector{Int},
    model::Vector{Int},
    wl::Int, ws::Int, dm::Int, order::Int, nr_tau::Int,
)::CTResult
    n_pairs = length(pairs)
    n_win = length(pairs[1].timepoints)
    n_coeff = isempty(pairs[1].timepoints) ? 0 : length(pairs[1].timepoints[1].coefficients)

    coefficients = Array{Float64,3}(undef, n_pairs, n_win, n_coeff)
    errors = Matrix{Float64}(undef, n_pairs, n_win)
    win_starts = Vector{Int64}(undef, n_win)
    win_ends = Vector{Int64}(undef, n_win)

    for (pi, pair_data) in enumerate(pairs)
        for (wi, tp) in enumerate(pair_data.timepoints)
            if pi == 1
                win_starts[wi] = tp.window_start
                win_ends[wi] = tp.window_end
            end
            for (ki, c) in enumerate(tp.coefficients)
                coefficients[pi, wi, ki] = c
            end
            errors[pi, wi] = tp.error
        end
    end

    params = _make_params(;
        sfreq=sfreq, delays=delays, model=model,
        wl=wl, ws=ws, dm=dm, order=order, nr_tau=nr_tau,
    )
    return CTResult(coefficients, errors, win_starts, win_ends, pair_labels, params)
end

# =============================================================================
# run_de
# =============================================================================

"""
    run_de(data; sfreq=1.0, delays=..., kwargs...) -> DEResult

Run Dynamical Ergodicity analysis.

# Returns
A [`DEResult`](@ref) with an ergodicity measure per time window.
"""
function run_de(
    data::Matrix{Float64};
    sfreq::Float64=1.0,
    delays::Vector{Int}=collect(DDADefaults.DELAYS),
    model::Vector{Int}=copy(DDADefaults.MODEL_PARAMS),
    wl::Int=DDADefaults.WINDOW_LENGTH,
    ws::Int=DDADefaults.WINDOW_STEP,
    ct_wl::Union{Int,Nothing}=nothing,
    ct_ws::Union{Int,Nothing}=nothing,
    channel_labels::Union{Vector{String},Nothing}=nothing,
    binary_path::Union{String,Nothing}=nothing,
    dm::Int=DDADefaults.MODEL_DIMENSION,
    order::Int=DDADefaults.POLYNOMIAL_ORDER,
    nr_tau::Int=DDADefaults.NUM_TAU,
)::DEResult
    n_ch, _ = size(data)

    tmp = _write_temp_ascii(data)
    try
        channels_0based = collect(0:(n_ch - 1))
        request = DDARequest(
            tmp, channels_0based, ["DE"];
            window_length=wl, window_step=ws,
            ct_window_length=something(ct_wl, wl),
            ct_window_step=something(ct_ws, ws),
            delays=delays, model_encoding=model,
            dm=dm, order=order, nr_tau=nr_tau,
            sampling_rate=sfreq > 1000.0 ? sfreq : nothing,
        )
        runner = binary_path === nothing ? DDARunner() : DDARunner(binary_path)
        raw = run_analysis_structured(runner, request)

        haskey(raw, "DE") || error("No DE results in DDA output")
        return _raw_to_de(raw["DE"], sfreq, delays, model, wl, ws, dm, order, nr_tau)
    finally
        isfile(tmp) && rm(tmp; force=true)
    end
end

function _raw_to_de(
    channels::Vector{StructuredChannelData},
    sfreq::Float64,
    delays::Vector{Int},
    model::Vector{Int},
    wl::Int, ws::Int, dm::Int, order::Int, nr_tau::Int,
)::DEResult
    # DE produces a single value per window (stride=1, no coefficients)
    ch_data = channels[1]
    n_win = length(ch_data.timepoints)

    ergodicity = Vector{Float64}(undef, n_win)
    win_starts = Vector{Int64}(undef, n_win)
    win_ends = Vector{Int64}(undef, n_win)

    for (wi, tp) in enumerate(ch_data.timepoints)
        win_starts[wi] = tp.window_start
        win_ends[wi] = tp.window_end
        ergodicity[wi] = tp.error  # DE: the single value is stored as "error"
    end

    params = _make_params(;
        sfreq=sfreq, delays=delays, model=model,
        wl=wl, ws=ws, dm=dm, order=order, nr_tau=nr_tau,
    )
    return DEResult(ergodicity, win_starts, win_ends, params)
end

end # module API
