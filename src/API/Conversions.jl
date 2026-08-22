function _resolve_model_configuration(
    model::Runner.OptionalModelSpec,
    derivative_points::Union{Int, Nothing},
    dm::Union{Int, Nothing},
    order::Union{Int, Nothing},
    nr_tau::Int,
)::Tuple{Vector{Int}, Int, Int}
    Runner._validate_custom_model_request(
        model,
        nothing,
        derivative_points,
        dm,
        order,
    )
    resolved_derivative_points = Runner._resolve_derivative_points(
        derivative_points,
        dm,
    )
    resolved_order = something(order, DDADefaults.POLYNOMIAL_ORDER)
    resolved_model = Runner._resolve_model_terms(
        model,
        nothing,
        nr_tau,
        resolved_order,
    )
    return resolved_model, resolved_derivative_points, resolved_order
end

function _make_params(;
    delays::AbstractVector{<:Integer},
    model::AbstractVector{<:Integer},
    WL::Union{Int, Nothing},
    WS::Union{Int, Nothing},
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
        "delays" => Int[delays...],
        "model" => Int[model...],
        "WL" => WL,
        "WS" => WS,
        "derivative_points" => derivative_points,
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
    WL::Union{Int, Nothing},
    WS::Union{Int, Nothing};
    channels::Vector{StructuredChannelData}=StructuredChannelData[],
    start_offset::Int=0,
)::Tuple{Vector{Int64}, Vector{Int64}}
    if WL === nothing || WS === nothing
        return Runner._raw_window_bounds(channels)
    end
    return Runner._fixed_window_bounds(n_windows, WL, WS, start_offset)
end

function _time_axes(
    channels::Vector{StructuredChannelData},
    derivative_points::Int,
    TM::Int,
    sampling_rate,
    WL::Union{Int, Nothing},
    WS::Union{Int, Nothing};
    start_offset::Int=0,
)::Tuple{Vector{Float64}, Vector{Float64}, Vector{Int64}, Vector{Int64}}
    n_windows = isempty(channels) ? 0 : length(channels[1].timepoints)
    T = Runner._extract_raw_T(channels)
    t = Runner._compute_t_axis(T, derivative_points, TM, sampling_rate)
    window_starts, window_ends = _window_bounds(
        n_windows,
        WL,
        WS;
        channels=channels,
        start_offset=start_offset,
    )
    return T, t, window_starts, window_ends
end

function _st_from_raw(
    channels::Vector{StructuredChannelData},
    labels::Vector{String};
    delays::AbstractVector{<:Integer},
    model::AbstractVector{<:Integer},
    WL::Union{Int, Nothing},
    WS::Union{Int, Nothing},
    derivative_points::Int,
    TM::Int,
    order::Int,
    nr_tau::Int,
    sampling_rate,
    out_fn::Union{String, Nothing},
    selected_channels::Vector{Int},
)::STResult
    coefficients, errors = Runner._coefficient_arrays(channels)
    T, t, win_starts, win_ends = _time_axes(
        channels,
        derivative_points,
        TM,
        sampling_rate,
        WL,
        WS,
    )

    params = _make_params(;
        delays=delays,
        model=model,
        WL=WL,
        WS=WS,
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
    delays::AbstractVector{<:Integer},
    model::AbstractVector{<:Integer},
    WL::Union{Int, Nothing},
    WS::Union{Int, Nothing},
    derivative_points::Int,
    TM::Int,
    order::Int,
    nr_tau::Int,
    sampling_rate,
    out_fn::Union{String, Nothing},
    selected_channels::Vector{Int},
)::CTResult
    coefficients, errors = Runner._coefficient_arrays(pairs)
    T, t, win_starts, win_ends = _time_axes(
        pairs,
        derivative_points,
        TM,
        sampling_rate,
        WL,
        WS,
    )

    params = _make_params(;
        delays=delays,
        model=model,
        WL=WL,
        WS=WS,
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
    delays::AbstractVector{<:Integer},
    model::AbstractVector{<:Integer},
    WL::Union{Int, Nothing},
    WS::Union{Int, Nothing},
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
        WL,
        WS,
    )

    for (wi, tp) in enumerate(ch_data.timepoints)
        ergodicity[wi] = tp.value
    end

    params = _make_params(;
        delays=delays,
        model=model,
        WL=WL,
        WS=WS,
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
