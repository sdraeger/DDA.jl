# =============================================================================
# FILE-BASED API
# =============================================================================

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
    selected_channels = Runner._normalize_channels(channels)
    labels = _resolve_labels(file_path, selected_channels, channel_labels)
    model_terms, derivative_points_value, order_value = _resolve_model_configuration(
        model,
        derivative_points,
        dm,
        order,
        nr_tau,
    )
    sampling_rate_value = Runner._normalize_sampling_rate(sampling_rate)
    tm_value = Runner._resolve_tm(Int[delays...], TM)

    raw = run_analysis_structured(;
        file_path=file_path,
        channels=selected_channels,
        flavors=["ST"],
        WL=WL,
        WS=WS,
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
        WL=WL,
        WS=WS,
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
    selected_channels = Runner._normalize_channels(channels)
    length(selected_channels) >= 2 || error("CT analysis requires at least 2 channels, got $(length(selected_channels))")

    labels = _resolve_labels(file_path, selected_channels, channel_labels)
    pair_labels = _pair_labels(labels)
    model_terms, derivative_points_value, order_value = _resolve_model_configuration(
        model,
        derivative_points,
        dm,
        order,
        nr_tau,
    )
    sampling_rate_value = Runner._normalize_sampling_rate(sampling_rate)
    tm_value = Runner._resolve_tm(Int[delays...], TM)
    raw_pairs = StructuredChannelData[]

    for pair_channels in _pair_channel_sets(selected_channels)
        raw = run_analysis_structured(;
            file_path=file_path,
            channels=pair_channels,
            flavors=["CT"],
            WL=WL,
            WS=WS,
            ct_window_length=ct_wl,
            ct_window_step=ct_ws,
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
        WL=ct_wl,
        WS=ct_ws,
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
    selected_channels = Runner._normalize_channels(channels)
    model_terms, derivative_points_value, order_value = _resolve_model_configuration(
        model,
        derivative_points,
        dm,
        order,
        nr_tau,
    )
    sampling_rate_value = Runner._normalize_sampling_rate(sampling_rate)
    tm_value = Runner._resolve_tm(Int[delays...], TM)

    raw = run_analysis_structured(;
        file_path=file_path,
        channels=selected_channels,
        flavors=["DE"],
        WL=WL,
        WS=WS,
        ct_window_length=ct_wl,
        ct_window_step=ct_ws,
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
        WL=ct_wl,
        WS=ct_ws,
        derivative_points=derivative_points_value,
        TM=tm_value,
        order=order_value,
        nr_tau=nr_tau,
        sampling_rate=sampling_rate_value,
        out_fn=out_fn,
        selected_channels=selected_channels,
    )
end
