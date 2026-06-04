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
    sfreq::Float64,
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
        "sfreq" => sfreq,
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
        if isempty(channels)
            return (Int64[], Int64[])
        end
        starts = Int64[round(Int64, tp.window_start) for tp in channels[1].timepoints]
        ends = Int64[round(Int64, tp.window_end) for tp in channels[1].timepoints]
        return (starts, ends)
    end

    window_starts = Vector{Int64}(undef, n_windows)
    window_ends = Vector{Int64}(undef, n_windows)

    for window_idx in 1:n_windows
        window_start = Int64(start_offset + (window_idx - 1) * WS)
        window_end = Int64(window_start + WL)
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

function _coefficient_arrays(
    channels::Vector{StructuredChannelData},
)::Tuple{Array{Float64, 3}, Matrix{Float64}}
    n_entities = length(channels)
    n_windows = length(channels[1].timepoints)
    n_coefficients = isempty(channels[1].timepoints) ? 0 : length(channels[1].timepoints[1].coefficients)
    coefficients = Array{Float64, 3}(undef, n_entities, n_windows, n_coefficients)
    errors = Matrix{Float64}(undef, n_entities, n_windows)

    for (entity_idx, channel_data) in enumerate(channels)
        for (window_idx, tp) in enumerate(channel_data.timepoints)
            for (coefficient_idx, coeff) in enumerate(tp.coefficients)
                coefficients[entity_idx, window_idx, coefficient_idx] = coeff
            end
            errors[entity_idx, window_idx] = tp.error
        end
    end

    return coefficients, errors
end

function _st_from_raw(
    channels::Vector{StructuredChannelData},
    labels::Vector{String};
    sfreq::Float64,
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
    coefficients, errors = _coefficient_arrays(channels)
    T, t, win_starts, win_ends = _time_axes(
        channels,
        derivative_points,
        TM,
        sampling_rate,
        WL,
        WS,
    )

    params = _make_params(;
        sfreq=sfreq,
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
    sfreq::Float64,
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
    coefficients, errors = _coefficient_arrays(pairs)
    T, t, win_starts, win_ends = _time_axes(
        pairs,
        derivative_points,
        TM,
        sampling_rate,
        WL,
        WS,
    )

    params = _make_params(;
        sfreq=sfreq,
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
    sfreq::Float64,
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
        ergodicity[wi] = tp.error
    end

    params = _make_params(;
        sfreq=sfreq,
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
