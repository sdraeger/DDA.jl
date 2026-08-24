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
    variant::VariantResultData,
    WL::Union{Int, Nothing},
    WS::Union{Int, Nothing},
)::Tuple{Vector{Int64}, Vector{Int64}}
    if WL === nothing || WS === nothing
        return (variant.window_starts, variant.window_ends)
    end
    return Runner._fixed_window_bounds(size(variant.coefficients, 2), WL, WS)
end

function _time_axes(
    variant::VariantResultData,
    WL::Union{Int, Nothing},
    WS::Union{Int, Nothing},
)::Tuple{Vector{Float64}, Vector{Float64}, Vector{Int64}, Vector{Int64}}
    T = Float64.(variant.T[:, 1])
    window_starts, window_ends = _window_bounds(variant, WL, WS)
    return T, variant.t, window_starts, window_ends
end

function _st_from_raw(
    variant::VariantResultData,
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
    T, t, win_starts, win_ends = _time_axes(variant, WL, WS)

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
    return STResult(variant.coefficients, variant.errors, T, t, win_starts, win_ends, labels, params)
end

function _ct_from_raw(
    variants::Vector{VariantResultData},
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
    isempty(variants) && error("No CT results found in DDA output")
    coefficients = cat((v.coefficients for v in variants)...; dims=1)
    errors = vcat((v.errors for v in variants)...)
    T, t, win_starts, win_ends = _time_axes(first(variants), WL, WS)

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
    variant::VariantResultData;
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
    ergodicity = vec(variant.errors)
    T, t, win_starts, win_ends = _time_axes(variant, WL, WS)

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
