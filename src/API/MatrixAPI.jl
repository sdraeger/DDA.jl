# =============================================================================
# MATRIX-BASED API
# =============================================================================

"""
    run_st(; data, kwargs...) -> STResult

Run single-timeseries DDA on an in-memory channels × samples matrix.
"""
function _run_st_matrix(
    data::AbstractMatrix{<:Real};
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
    n_ch, _ = size(data)
    tmp = _write_temp_ascii(data)
    try
        return _run_st_file(
            tmp,
            collect(1:n_ch);
            sfreq=sfreq,
            delays=delays,
            model=model,
            WL=WL,
            WS=WS,
            channel_labels=channel_labels,
            binary_path=binary_path,
            out_fn=out_fn,
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
    run_ct(; data, kwargs...) -> CTResult

Run cross-timeseries DDA on an in-memory channels × samples matrix.
"""
function _run_ct_matrix(
    data::AbstractMatrix{<:Real};
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
    n_ch, _ = size(data)
    tmp = _write_temp_ascii(data)
    try
        return _run_ct_file(
            tmp,
            collect(1:n_ch);
            sfreq=sfreq,
            delays=delays,
            model=model,
            WL=WL,
            WS=WS,
            ct_wl=ct_wl,
            ct_ws=ct_ws,
            channel_labels=channel_labels,
            binary_path=binary_path,
            out_fn=out_fn,
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
    run_de(; data, kwargs...) -> DEResult

Run dynamical-ergodicity DDA on an in-memory channels × samples matrix.
"""
function _run_de_matrix(
    data::AbstractMatrix{<:Real};
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
    n_ch, _ = size(data)
    tmp = _write_temp_ascii(data)
    try
        return _run_de_file(
            tmp,
            collect(1:n_ch);
            sfreq=sfreq,
            delays=delays,
            model=model,
            WL=WL,
            WS=WS,
            ct_wl=ct_wl,
            ct_ws=ct_ws,
            binary_path=binary_path,
            out_fn=out_fn,
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
