# =============================================================================
# FILE-BASED API
# =============================================================================

function _file_run_context(
    file_path::AbstractString,
    channels::AbstractVector{<:Integer};
    channel_labels::Union{Vector{String}, Nothing},
    model::Runner.OptionalModelSpec,
    derivative_points::Union{Int, Nothing},
    dm::Union{Int, Nothing},
    order::Union{Int, Nothing},
    nr_tau::Int,
    sampling_rate,
    delays::AbstractVector{<:Integer},
    TM::Union{Int, Nothing},
    resolve_labels::Bool=true,
)
    selected_channels = Runner._normalize_channels(channels)
    labels = resolve_labels ? _resolve_labels(file_path, selected_channels, channel_labels) : String[]
    model_terms, derivative_points_value, order_value = _resolve_model_configuration(
        model,
        derivative_points,
        dm,
        order,
        nr_tau,
    )
    return (
        selected_channels=selected_channels,
        labels=labels,
        model_terms=model_terms,
        derivative_points=derivative_points_value,
        order=order_value,
        sampling_rate=Runner._normalize_sampling_rate(sampling_rate),
        TM=Runner._resolve_tm(Int[delays...], TM),
    )
end

function _run_file_flavor(
    file_path::AbstractString,
    selected_channels::Vector{Int},
    flavor::String;
    WL::Union{Int, Nothing},
    WS::Union{Int, Nothing},
    ct_wl::Union{Int, Nothing}=nothing,
    ct_ws::Union{Int, Nothing}=nothing,
    delays::AbstractVector{<:Integer},
    model_terms::AbstractVector{<:Integer},
    derivative_points::Int,
    order::Int,
    nr_tau::Int,
    sampling_rate,
    TM::Int,
    out_fn::Union{String, Nothing},
    binary_path::Union{String, Nothing},
)
    return run_analysis_structured(;
        file_path=file_path,
        channels=selected_channels,
        flavors=[flavor],
        WL=WL,
        WS=WS,
        ct_window_length=ct_wl,
        ct_window_step=ct_ws,
        delays=Int[delays...],
        model=model_terms,
        derivative_points=derivative_points,
        order=order,
        nr_tau=nr_tau,
        sampling_rate=sampling_rate,
        TM=TM,
        out_fn=out_fn,
        binary_path=binary_path,
    )
end

"""
    run_st(; file_path, channels, kwargs...) -> STResult

Run single-timeseries DDA directly on an EDF or ASCII file using 1-based channel indices.

Important keywords:
- `binary_path`: resolve the DDA binary without relying on environment variables
- `model`: optional custom model passed to `-MODEL` as indices, or as matrix rows mapped to indices
- `derivative_points`: preferred name for binary `-dm`
- Passing a custom `model` also requires explicit `derivative_points` and `order`
- `TM`: optional offset used only to compute `result.t`; defaults to `max(delays)`
- `out_fn`: optional `-OUT_FN` base; defaults to a temporary path for the call
- `sampling_rate`: optional value passed to `-SR`. A scalar emits `-SR N`; a tuple emits `-SR N1 N2`; `nothing` emits no `-SR`
"""
function _run_st_file(
    file_path::AbstractString,
    channels::AbstractVector{<:Integer};
    sfreq::Float64=1.0,
    delays::AbstractVector{<:Integer}=collect(DDADefaults.DELAYS),
    model::Runner.OptionalModelSpec=nothing,
    WL::Union{Int, Nothing}=DDADefaults.WL,
    WS::Union{Int, Nothing}=DDADefaults.WS,
    channel_labels::Union{Vector{String}, Nothing}=nothing,
    binary_path::Union{String, Nothing}=nothing,
    out_fn::Union{String, Nothing}=nothing,
    derivative_points::Union{Int, Nothing}=nothing,
    dm::Union{Int, Nothing}=nothing,
    order::Union{Int, Nothing}=nothing,
    nr_tau::Int=DDADefaults.NUM_TAU,
    sampling_rate=DDADefaults.SAMPLING_RATE,
    TM::Union{Int, Nothing}=nothing,
)::STResult
    ctx = _file_run_context(
        file_path,
        channels;
        channel_labels=channel_labels,
        model=model,
        derivative_points=derivative_points,
        dm=dm,
        order=order,
        nr_tau=nr_tau,
        sampling_rate=sampling_rate,
        delays=delays,
        TM=TM,
    )

    raw = _run_file_flavor(
        file_path,
        ctx.selected_channels,
        "ST";
        WL=WL,
        WS=WS,
        delays=delays,
        model_terms=ctx.model_terms,
        derivative_points=ctx.derivative_points,
        order=ctx.order,
        nr_tau=nr_tau,
        sampling_rate=ctx.sampling_rate,
        TM=ctx.TM,
        out_fn=out_fn,
        binary_path=binary_path,
    )

    haskey(raw, "ST") || error("No ST results found in DDA output")
    return _st_from_raw(
        raw["ST"],
        ctx.labels;
        sfreq=sfreq,
        delays=delays,
        model=ctx.model_terms,
        WL=WL,
        WS=WS,
        derivative_points=ctx.derivative_points,
        TM=ctx.TM,
        order=ctx.order,
        nr_tau=nr_tau,
        sampling_rate=ctx.sampling_rate,
        out_fn=out_fn,
        selected_channels=ctx.selected_channels,
    )
end

"""
    run_ct(; file_path, channels, kwargs...) -> CTResult

Run cross-timeseries DDA directly on an EDF or ASCII file using 1-based channel indices.
"""
function _run_ct_file(
    file_path::AbstractString,
    channels::AbstractVector{<:Integer};
    sfreq::Float64=1.0,
    delays::AbstractVector{<:Integer}=collect(DDADefaults.DELAYS),
    model::Runner.OptionalModelSpec=nothing,
    WL::Union{Int, Nothing}=DDADefaults.WL,
    WS::Union{Int, Nothing}=DDADefaults.WS,
    ct_wl::Union{Int, Nothing}=nothing,
    ct_ws::Union{Int, Nothing}=nothing,
    channel_labels::Union{Vector{String}, Nothing}=nothing,
    binary_path::Union{String, Nothing}=nothing,
    out_fn::Union{String, Nothing}=nothing,
    derivative_points::Union{Int, Nothing}=nothing,
    dm::Union{Int, Nothing}=nothing,
    order::Union{Int, Nothing}=nothing,
    nr_tau::Int=DDADefaults.NUM_TAU,
    sampling_rate=DDADefaults.SAMPLING_RATE,
    TM::Union{Int, Nothing}=nothing,
)::CTResult
    ctx = _file_run_context(
        file_path,
        channels;
        channel_labels=channel_labels,
        model=model,
        derivative_points=derivative_points,
        dm=dm,
        order=order,
        nr_tau=nr_tau,
        sampling_rate=sampling_rate,
        delays=delays,
        TM=TM,
    )
    length(ctx.selected_channels) >= 2 || error("CT analysis requires at least 2 channels, got $(length(ctx.selected_channels))")

    pair_labels = _pair_labels(ctx.labels)
    raw_pairs = StructuredChannelData[]

    for pair_channels in _pair_channel_sets(ctx.selected_channels)
        raw = _run_file_flavor(
            file_path,
            pair_channels,
            "CT";
            WL=WL,
            WS=WS,
            ct_wl=ct_wl,
            ct_ws=ct_ws,
            delays=delays,
            model_terms=ctx.model_terms,
            derivative_points=ctx.derivative_points,
            order=ctx.order,
            nr_tau=nr_tau,
            sampling_rate=ctx.sampling_rate,
            TM=ctx.TM,
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
        model=ctx.model_terms,
        WL=ct_wl,
        WS=ct_ws,
        derivative_points=ctx.derivative_points,
        TM=ctx.TM,
        order=ctx.order,
        nr_tau=nr_tau,
        sampling_rate=ctx.sampling_rate,
        out_fn=out_fn,
        selected_channels=ctx.selected_channels,
    )
end

"""
    run_de(; file_path, channels, kwargs...) -> DEResult

Run dynamical-ergodicity DDA directly on an EDF or ASCII file using 1-based channel indices.
"""
function _run_de_file(
    file_path::AbstractString,
    channels::AbstractVector{<:Integer};
    sfreq::Float64=1.0,
    delays::AbstractVector{<:Integer}=collect(DDADefaults.DELAYS),
    model::Runner.OptionalModelSpec=nothing,
    WL::Union{Int, Nothing}=DDADefaults.WL,
    WS::Union{Int, Nothing}=DDADefaults.WS,
    ct_wl::Union{Int, Nothing}=nothing,
    ct_ws::Union{Int, Nothing}=nothing,
    binary_path::Union{String, Nothing}=nothing,
    out_fn::Union{String, Nothing}=nothing,
    derivative_points::Union{Int, Nothing}=nothing,
    dm::Union{Int, Nothing}=nothing,
    order::Union{Int, Nothing}=nothing,
    nr_tau::Int=DDADefaults.NUM_TAU,
    sampling_rate=DDADefaults.SAMPLING_RATE,
    TM::Union{Int, Nothing}=nothing,
)::DEResult
    ctx = _file_run_context(
        file_path,
        channels;
        channel_labels=nothing,
        model=model,
        derivative_points=derivative_points,
        dm=dm,
        order=order,
        nr_tau=nr_tau,
        sampling_rate=sampling_rate,
        delays=delays,
        TM=TM,
        resolve_labels=false,
    )

    raw = _run_file_flavor(
        file_path,
        ctx.selected_channels,
        "DE";
        WL=WL,
        WS=WS,
        ct_wl=ct_wl,
        ct_ws=ct_ws,
        delays=delays,
        model_terms=ctx.model_terms,
        derivative_points=ctx.derivative_points,
        order=ctx.order,
        nr_tau=nr_tau,
        sampling_rate=ctx.sampling_rate,
        TM=ctx.TM,
        out_fn=out_fn,
        binary_path=binary_path,
    )

    haskey(raw, "DE") || error("No DE results found in DDA output")
    return _de_from_raw(
        raw["DE"];
        sfreq=sfreq,
        delays=delays,
        model=ctx.model_terms,
        WL=ct_wl,
        WS=ct_ws,
        derivative_points=ctx.derivative_points,
        TM=ctx.TM,
        order=ctx.order,
        nr_tau=nr_tau,
        sampling_rate=ctx.sampling_rate,
        out_fn=out_fn,
        selected_channels=ctx.selected_channels,
    )
end
