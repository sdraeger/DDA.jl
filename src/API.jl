"""High-level analysis functions for the DDA binary."""
module API

using Printf
using ..DDADefaults
using ..Results
using ..Runner

export run_st, run_ct, run_de

# =============================================================================
# HELPERS
# =============================================================================

"""Write a channels × samples matrix to a temporary ASCII file for the DDA binary."""
function _write_temp_ascii(data::AbstractMatrix{<:Real})::String
    path = tempname() * ".txt"
    open(path, "w") do io
        n_ch, n_samp = size(data)
        for t in 1:n_samp
            for ch in 1:n_ch
                ch > 1 && print(io, '\t')
                @printf(io, "%.15f", Float64(data[ch, t]))
            end
            println(io)
        end
    end
    return path
end

function _resolve_labels(
    channels::AbstractVector{<:Integer},
    channel_labels::Union{Vector{String}, Nothing},
)::Vector{String}
    labels = channel_labels === nothing ? ["ch$(ch)" for ch in channels] : copy(channel_labels)
    length(labels) == length(channels) || error("Expected $(length(channels)) channel labels, got $(length(labels))")
    return labels
end

function _resolve_labels(
    file_path::AbstractString,
    channels::AbstractVector{<:Integer},
    channel_labels::Union{Vector{String}, Nothing},
)::Vector{String}
    labels = if channel_labels === nothing
        Runner._resolve_requested_channel_labels(file_path, channels; fallback_prefix="ch")
    else
        copy(channel_labels)
    end
    length(labels) == length(channels) || error("Expected $(length(channels)) channel labels, got $(length(labels))")
    return labels
end

function _pair_labels(labels::Vector{String})::Vector{String}
    pairs = String[]
    for i in 1:length(labels), j in (i + 1):length(labels)
        push!(pairs, "$(labels[i])-$(labels[j])")
    end
    return pairs
end

function _pair_channel_sets(channels::Vector{Int})::Vector{Vector{Int}}
    pairs = Vector{Vector{Int}}()
    for i in 1:length(channels), j in (i + 1):length(channels)
        push!(pairs, [channels[i], channels[j]])
    end
    return pairs
end

function _resolve_model_configuration(
    model::Union{AbstractVector{<:Integer}, Nothing},
    model_dimension::Union{Int, Nothing},
    derivative_points::Union{Int, Nothing},
    dm::Union{Int, Nothing},
    order::Union{Int, Nothing},
)::Tuple{Vector{Int}, Int, Int}
    Runner._validate_custom_model_request(
        model,
        nothing,
        model_dimension,
        derivative_points,
        dm,
        order,
    )
    resolved_model = Int[something(model, copy(DDADefaults.MODEL_PARAMS))...]
    resolved_derivative_points = Runner._resolve_derivative_points(
        model_dimension,
        derivative_points,
        dm,
    )
    resolved_order = something(order, DDADefaults.POLYNOMIAL_ORDER)
    return resolved_model, resolved_derivative_points, resolved_order
end

function _make_params(;
    sfreq::Float64,
    delays::AbstractVector{<:Integer},
    model::AbstractVector{<:Integer},
    wl::Int,
    ws::Int,
    derivative_points::Int,
    TM::Int,
    order::Int,
    nr_tau::Int,
    sampling_rate,
    channels::Union{Vector{Int}, Nothing}=nothing,
    out_fn::Union{String, Nothing}=nothing,
    extra::Dict{String, Any}=Dict{String, Any}(),
)::Dict{String, Any}
    params = Dict{String, Any}(
        "sfreq" => sfreq,
        "delays" => Int[delays...],
        "model" => Int[model...],
        "wl" => wl,
        "ws" => ws,
        "derivative_points" => derivative_points,
        "model_dimension" => derivative_points,
        "dm" => derivative_points,
        "TM" => TM,
        "order" => order,
        "nr_tau" => nr_tau,
        "sampling_rate" => sampling_rate,
        "out_fn" => out_fn,
    )
    channels !== nothing && (params["channels"] = copy(channels))
    merge!(params, extra)
    return params
end

function _window_bounds(
    n_windows::Int,
    window_length::Int,
    window_step::Int;
    start_offset::Int=0,
)::Tuple{Vector{Int64}, Vector{Int64}}
    window_starts = Vector{Int64}(undef, n_windows)
    window_ends = Vector{Int64}(undef, n_windows)

    for window_idx in 1:n_windows
        window_start = Int64(start_offset + (window_idx - 1) * window_step)
        window_end = Int64(window_start + window_length)
        window_starts[window_idx] = window_start
        window_ends[window_idx] = window_end
    end

    return window_starts, window_ends
end

function _time_axes(
    channels::Vector{StructuredChannelData},
    derivative_points::Int,
    TM::Int,
    sampling_rate,
    window_length::Int,
    window_step::Int;
    start_offset::Int=0,
)::Tuple{Vector{Float64}, Vector{Float64}, Vector{Int64}, Vector{Int64}}
    n_windows = isempty(channels) ? 0 : length(channels[1].timepoints)
    T = Runner._extract_raw_T(channels)
    t = Runner._compute_t_axis(T, derivative_points, TM, sampling_rate)
    window_starts, window_ends = _window_bounds(
        n_windows,
        window_length,
        window_step;
        start_offset=start_offset,
    )
    return T, t, window_starts, window_ends
end

function _st_from_raw(
    channels::Vector{StructuredChannelData},
    labels::Vector{String};
    sfreq::Float64,
    delays::AbstractVector{<:Integer},
    model::AbstractVector{<:Integer},
    wl::Int,
    ws::Int,
    derivative_points::Int,
    TM::Int,
    order::Int,
    nr_tau::Int,
    sampling_rate,
    out_fn::Union{String, Nothing},
    selected_channels::Vector{Int},
)::STResult
    n_ch = length(channels)
    n_win = length(channels[1].timepoints)
    n_coeff = isempty(channels[1].timepoints) ? 0 : length(channels[1].timepoints[1].coefficients)

    coefficients = Array{Float64, 3}(undef, n_ch, n_win, n_coeff)
    errors = Matrix{Float64}(undef, n_ch, n_win)
    T, t, win_starts, win_ends = _time_axes(
        channels,
        derivative_points,
        TM,
        sampling_rate,
        wl,
        ws,
    )

    for (ci, ch_data) in enumerate(channels)
        for (wi, tp) in enumerate(ch_data.timepoints)
            for (ki, coeff) in enumerate(tp.coefficients)
                coefficients[ci, wi, ki] = coeff
            end
            errors[ci, wi] = tp.error
        end
    end

    params = _make_params(;
        sfreq=sfreq,
        delays=delays,
        model=model,
        wl=wl,
        ws=ws,
        derivative_points=derivative_points,
        TM=TM,
        order=order,
        nr_tau=nr_tau,
        sampling_rate=sampling_rate,
        channels=selected_channels,
        out_fn=out_fn,
    )
    return STResult(coefficients, errors, T, t, win_starts, win_ends, labels, params)
end

function _ct_from_raw(
    pairs::Vector{StructuredChannelData},
    pair_labels::Vector{String};
    sfreq::Float64,
    delays::AbstractVector{<:Integer},
    model::AbstractVector{<:Integer},
    wl::Int,
    ws::Int,
    derivative_points::Int,
    TM::Int,
    order::Int,
    nr_tau::Int,
    sampling_rate,
    out_fn::Union{String, Nothing},
    selected_channels::Vector{Int},
)::CTResult
    n_pairs = length(pairs)
    n_win = length(pairs[1].timepoints)
    n_coeff = isempty(pairs[1].timepoints) ? 0 : length(pairs[1].timepoints[1].coefficients)

    coefficients = Array{Float64, 3}(undef, n_pairs, n_win, n_coeff)
    errors = Matrix{Float64}(undef, n_pairs, n_win)
    T, t, win_starts, win_ends = _time_axes(
        pairs,
        derivative_points,
        TM,
        sampling_rate,
        wl,
        ws,
    )

    for (pi, pair_data) in enumerate(pairs)
        for (wi, tp) in enumerate(pair_data.timepoints)
            for (ki, coeff) in enumerate(tp.coefficients)
                coefficients[pi, wi, ki] = coeff
            end
            errors[pi, wi] = tp.error
        end
    end

    params = _make_params(;
        sfreq=sfreq,
        delays=delays,
        model=model,
        wl=wl,
        ws=ws,
        derivative_points=derivative_points,
        TM=TM,
        order=order,
        nr_tau=nr_tau,
        sampling_rate=sampling_rate,
        channels=selected_channels,
        out_fn=out_fn,
    )
    return CTResult(coefficients, errors, T, t, win_starts, win_ends, pair_labels, params)
end

function _de_from_raw(
    channels::Vector{StructuredChannelData};
    sfreq::Float64,
    delays::AbstractVector{<:Integer},
    model::AbstractVector{<:Integer},
    wl::Int,
    ws::Int,
    derivative_points::Int,
    TM::Int,
    order::Int,
    nr_tau::Int,
    sampling_rate,
    out_fn::Union{String, Nothing},
    selected_channels::Vector{Int},
)::DEResult
    ch_data = channels[1]
    n_win = length(ch_data.timepoints)

    ergodicity = Vector{Float64}(undef, n_win)
    T, t, win_starts, win_ends = _time_axes(
        channels,
        derivative_points,
        TM,
        sampling_rate,
        wl,
        ws,
    )

    for (wi, tp) in enumerate(ch_data.timepoints)
        ergodicity[wi] = tp.error
    end

    params = _make_params(;
        sfreq=sfreq,
        delays=delays,
        model=model,
        wl=wl,
        ws=ws,
        derivative_points=derivative_points,
        TM=TM,
        order=order,
        nr_tau=nr_tau,
        sampling_rate=sampling_rate,
        channels=selected_channels,
        out_fn=out_fn,
    )
    return DEResult(ergodicity, T, t, win_starts, win_ends, params)
end

# =============================================================================
# FILE-BASED API
# =============================================================================

"""
    run_st(file_path, channels; kwargs...) -> STResult

Run single-timeseries DDA directly on an EDF or ASCII file using 1-based channel indices.

Important keywords:
- `binary_path`: resolve the DDA binary without relying on environment variables
- `model`: optional custom values passed to `-MODEL`
- `derivative_points`: preferred name for binary `-dm`; `model_dimension` remains a compatibility alias
- Passing a custom `model` also requires explicit `model_dimension` or `derivative_points`, and `order`
- `TM`: optional offset used only to compute `result.t`; defaults to `max(delays)`
- `out_fn`: optional `-OUT_FN` base; defaults to a temporary path for the call
- `sampling_rate`: optional values passed to `-SR`. If you pass `(N, N)`, it is only used for `result.t`
"""
function run_st(
    file_path::AbstractString,
    channels::AbstractVector{<:Integer};
    sfreq::Float64=1.0,
    delays::AbstractVector{<:Integer}=collect(DDADefaults.DELAYS),
    model::Union{Vector{Int}, Nothing}=nothing,
    wl::Int=DDADefaults.WINDOW_LENGTH,
    ws::Int=DDADefaults.WINDOW_STEP,
    channel_labels::Union{Vector{String}, Nothing}=nothing,
    binary_path::Union{String, Nothing}=nothing,
    out_fn::Union{String, Nothing}=nothing,
    model_dimension::Union{Int, Nothing}=nothing,
    derivative_points::Union{Int, Nothing}=nothing,
    dm::Union{Int, Nothing}=nothing,
    order::Union{Int, Nothing}=nothing,
    nr_tau::Int=DDADefaults.NUM_TAU,
    sampling_rate=DDADefaults.SAMPLING_RATE,
    TM::Union{Int, Nothing}=nothing,
)::STResult
    selected_channels = Runner._normalize_channels(channels)
    labels = _resolve_labels(file_path, selected_channels, channel_labels)
    model_terms, derivative_points_value, order_value = _resolve_model_configuration(
        model,
        model_dimension,
        derivative_points,
        dm,
        order,
    )
    sampling_rate_value = Runner._normalize_sampling_rate(sampling_rate)
    tm_value = Runner._resolve_tm(Int[delays...], TM)

    raw = run_analysis_structured(
        file_path,
        selected_channels,
        ["ST"];
        window_length=wl,
        window_step=ws,
        delays=Int[delays...],
        model=model_terms,
        derivative_points=derivative_points_value,
        order=order_value,
        nr_tau=nr_tau,
        sampling_rate=sampling_rate_value,
        TM=tm_value,
        out_fn=out_fn,
        binary_path=binary_path,
    )

    haskey(raw, "ST") || error("No ST results found in DDA output")
    return _st_from_raw(
        raw["ST"],
        labels;
        sfreq=sfreq,
        delays=delays,
        model=model_terms,
        wl=wl,
        ws=ws,
        derivative_points=derivative_points_value,
        TM=tm_value,
        order=order_value,
        nr_tau=nr_tau,
        sampling_rate=sampling_rate_value,
        out_fn=out_fn,
        selected_channels=selected_channels,
    )
end

"""
    run_ct(file_path, channels; kwargs...) -> CTResult

Run cross-timeseries DDA directly on an EDF or ASCII file using 1-based channel indices.
"""
function run_ct(
    file_path::AbstractString,
    channels::AbstractVector{<:Integer};
    sfreq::Float64=1.0,
    delays::AbstractVector{<:Integer}=collect(DDADefaults.DELAYS),
    model::Union{Vector{Int}, Nothing}=nothing,
    wl::Int=DDADefaults.WINDOW_LENGTH,
    ws::Int=DDADefaults.WINDOW_STEP,
    ct_wl::Union{Int, Nothing}=nothing,
    ct_ws::Union{Int, Nothing}=nothing,
    channel_labels::Union{Vector{String}, Nothing}=nothing,
    binary_path::Union{String, Nothing}=nothing,
    out_fn::Union{String, Nothing}=nothing,
    model_dimension::Union{Int, Nothing}=nothing,
    derivative_points::Union{Int, Nothing}=nothing,
    dm::Union{Int, Nothing}=nothing,
    order::Union{Int, Nothing}=nothing,
    nr_tau::Int=DDADefaults.NUM_TAU,
    sampling_rate=DDADefaults.SAMPLING_RATE,
    TM::Union{Int, Nothing}=nothing,
)::CTResult
    selected_channels = Runner._normalize_channels(channels)
    length(selected_channels) >= 2 || error("CT analysis requires at least 2 channels, got $(length(selected_channels))")

    labels = _resolve_labels(file_path, selected_channels, channel_labels)
    pair_labels = _pair_labels(labels)
    model_terms, derivative_points_value, order_value = _resolve_model_configuration(
        model,
        model_dimension,
        derivative_points,
        dm,
        order,
    )
    sampling_rate_value = Runner._normalize_sampling_rate(sampling_rate)
    tm_value = Runner._resolve_tm(Int[delays...], TM)
    raw_pairs = StructuredChannelData[]

    for pair_channels in _pair_channel_sets(selected_channels)
        raw = run_analysis_structured(
            file_path,
            pair_channels,
            ["CT"];
            window_length=wl,
            window_step=ws,
            ct_window_length=something(ct_wl, wl),
            ct_window_step=something(ct_ws, ws),
            delays=Int[delays...],
            model=model_terms,
            derivative_points=derivative_points_value,
            order=order_value,
            nr_tau=nr_tau,
            sampling_rate=sampling_rate_value,
            TM=tm_value,
            out_fn=out_fn,
            binary_path=binary_path,
        )

        haskey(raw, "CT") || error("No CT results found in DDA output for channels $(pair_channels)")
        isempty(raw["CT"]) && error("Empty CT result returned for channels $(pair_channels)")
        push!(raw_pairs, raw["CT"][1])
    end

    return _ct_from_raw(
        raw_pairs,
        pair_labels;
        sfreq=sfreq,
        delays=delays,
        model=model_terms,
        wl=something(ct_wl, wl),
        ws=something(ct_ws, ws),
        derivative_points=derivative_points_value,
        TM=tm_value,
        order=order_value,
        nr_tau=nr_tau,
        sampling_rate=sampling_rate_value,
        out_fn=out_fn,
        selected_channels=selected_channels,
    )
end

"""
    run_de(file_path, channels; kwargs...) -> DEResult

Run dynamical-ergodicity DDA directly on an EDF or ASCII file using 1-based channel indices.
"""
function run_de(
    file_path::AbstractString,
    channels::AbstractVector{<:Integer};
    sfreq::Float64=1.0,
    delays::AbstractVector{<:Integer}=collect(DDADefaults.DELAYS),
    model::Union{Vector{Int}, Nothing}=nothing,
    wl::Int=DDADefaults.WINDOW_LENGTH,
    ws::Int=DDADefaults.WINDOW_STEP,
    ct_wl::Union{Int, Nothing}=nothing,
    ct_ws::Union{Int, Nothing}=nothing,
    binary_path::Union{String, Nothing}=nothing,
    out_fn::Union{String, Nothing}=nothing,
    model_dimension::Union{Int, Nothing}=nothing,
    derivative_points::Union{Int, Nothing}=nothing,
    dm::Union{Int, Nothing}=nothing,
    order::Union{Int, Nothing}=nothing,
    nr_tau::Int=DDADefaults.NUM_TAU,
    sampling_rate=DDADefaults.SAMPLING_RATE,
    TM::Union{Int, Nothing}=nothing,
)::DEResult
    selected_channels = Runner._normalize_channels(channels)
    model_terms, derivative_points_value, order_value = _resolve_model_configuration(
        model,
        model_dimension,
        derivative_points,
        dm,
        order,
    )
    sampling_rate_value = Runner._normalize_sampling_rate(sampling_rate)
    tm_value = Runner._resolve_tm(Int[delays...], TM)

    raw = run_analysis_structured(
        file_path,
        selected_channels,
        ["DE"];
        window_length=wl,
        window_step=ws,
        ct_window_length=something(ct_wl, wl),
        ct_window_step=something(ct_ws, ws),
        delays=Int[delays...],
        model=model_terms,
        derivative_points=derivative_points_value,
        order=order_value,
        nr_tau=nr_tau,
        sampling_rate=sampling_rate_value,
        TM=tm_value,
        out_fn=out_fn,
        binary_path=binary_path,
    )

    haskey(raw, "DE") || error("No DE results found in DDA output")
    return _de_from_raw(
        raw["DE"];
        sfreq=sfreq,
        delays=delays,
        model=model_terms,
        wl=something(ct_wl, wl),
        ws=something(ct_ws, ws),
        derivative_points=derivative_points_value,
        TM=tm_value,
        order=order_value,
        nr_tau=nr_tau,
        sampling_rate=sampling_rate_value,
        out_fn=out_fn,
        selected_channels=selected_channels,
    )
end

# =============================================================================
# MATRIX-BASED API
# =============================================================================

"""
    run_st(data; kwargs...) -> STResult

Run single-timeseries DDA on an in-memory channels × samples matrix.
"""
function run_st(
    data::AbstractMatrix{<:Real};
    sfreq::Float64=1.0,
    delays::AbstractVector{<:Integer}=collect(DDADefaults.DELAYS),
    model::Union{Vector{Int}, Nothing}=nothing,
    wl::Int=DDADefaults.WINDOW_LENGTH,
    ws::Int=DDADefaults.WINDOW_STEP,
    channel_labels::Union{Vector{String}, Nothing}=nothing,
    binary_path::Union{String, Nothing}=nothing,
    out_fn::Union{String, Nothing}=nothing,
    model_dimension::Union{Int, Nothing}=nothing,
    derivative_points::Union{Int, Nothing}=nothing,
    dm::Union{Int, Nothing}=nothing,
    order::Union{Int, Nothing}=nothing,
    nr_tau::Int=DDADefaults.NUM_TAU,
    sampling_rate=DDADefaults.SAMPLING_RATE,
    TM::Union{Int, Nothing}=nothing,
)::STResult
    n_ch, _ = size(data)
    tmp = _write_temp_ascii(data)
    try
        return run_st(
            tmp,
            collect(1:n_ch);
            sfreq=sfreq,
            delays=delays,
            model=model,
            wl=wl,
            ws=ws,
            channel_labels=channel_labels,
            binary_path=binary_path,
            out_fn=out_fn,
            model_dimension=model_dimension,
            derivative_points=derivative_points,
            dm=dm,
            order=order,
            nr_tau=nr_tau,
            sampling_rate=sampling_rate,
            TM=TM,
        )
    finally
        isfile(tmp) && rm(tmp; force=true)
    end
end

"""
    run_ct(data; kwargs...) -> CTResult

Run cross-timeseries DDA on an in-memory channels × samples matrix.
"""
function run_ct(
    data::AbstractMatrix{<:Real};
    sfreq::Float64=1.0,
    delays::AbstractVector{<:Integer}=collect(DDADefaults.DELAYS),
    model::Union{Vector{Int}, Nothing}=nothing,
    wl::Int=DDADefaults.WINDOW_LENGTH,
    ws::Int=DDADefaults.WINDOW_STEP,
    ct_wl::Union{Int, Nothing}=nothing,
    ct_ws::Union{Int, Nothing}=nothing,
    channel_labels::Union{Vector{String}, Nothing}=nothing,
    binary_path::Union{String, Nothing}=nothing,
    out_fn::Union{String, Nothing}=nothing,
    model_dimension::Union{Int, Nothing}=nothing,
    derivative_points::Union{Int, Nothing}=nothing,
    dm::Union{Int, Nothing}=nothing,
    order::Union{Int, Nothing}=nothing,
    nr_tau::Int=DDADefaults.NUM_TAU,
    sampling_rate=DDADefaults.SAMPLING_RATE,
    TM::Union{Int, Nothing}=nothing,
)::CTResult
    n_ch, _ = size(data)
    tmp = _write_temp_ascii(data)
    try
        return run_ct(
            tmp,
            collect(1:n_ch);
            sfreq=sfreq,
            delays=delays,
            model=model,
            wl=wl,
            ws=ws,
            ct_wl=ct_wl,
            ct_ws=ct_ws,
            channel_labels=channel_labels,
            binary_path=binary_path,
            out_fn=out_fn,
            model_dimension=model_dimension,
            derivative_points=derivative_points,
            dm=dm,
            order=order,
            nr_tau=nr_tau,
            sampling_rate=sampling_rate,
            TM=TM,
        )
    finally
        isfile(tmp) && rm(tmp; force=true)
    end
end

"""
    run_de(data; kwargs...) -> DEResult

Run dynamical-ergodicity DDA on an in-memory channels × samples matrix.
"""
function run_de(
    data::AbstractMatrix{<:Real};
    sfreq::Float64=1.0,
    delays::AbstractVector{<:Integer}=collect(DDADefaults.DELAYS),
    model::Union{Vector{Int}, Nothing}=nothing,
    wl::Int=DDADefaults.WINDOW_LENGTH,
    ws::Int=DDADefaults.WINDOW_STEP,
    ct_wl::Union{Int, Nothing}=nothing,
    ct_ws::Union{Int, Nothing}=nothing,
    binary_path::Union{String, Nothing}=nothing,
    out_fn::Union{String, Nothing}=nothing,
    model_dimension::Union{Int, Nothing}=nothing,
    derivative_points::Union{Int, Nothing}=nothing,
    dm::Union{Int, Nothing}=nothing,
    order::Union{Int, Nothing}=nothing,
    nr_tau::Int=DDADefaults.NUM_TAU,
    sampling_rate=DDADefaults.SAMPLING_RATE,
    TM::Union{Int, Nothing}=nothing,
)::DEResult
    n_ch, _ = size(data)
    tmp = _write_temp_ascii(data)
    try
        return run_de(
            tmp,
            collect(1:n_ch);
            sfreq=sfreq,
            delays=delays,
            model=model,
            wl=wl,
            ws=ws,
            ct_wl=ct_wl,
            ct_ws=ct_ws,
            binary_path=binary_path,
            out_fn=out_fn,
            model_dimension=model_dimension,
            derivative_points=derivative_points,
            dm=dm,
            order=order,
            nr_tau=nr_tau,
            sampling_rate=sampling_rate,
            TM=TM,
        )
    finally
        isfile(tmp) && rm(tmp; force=true)
    end
end

end # module API
