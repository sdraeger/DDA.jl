function _normalize_channels(channels::AbstractVector{<:Integer})::Vector{Int}
    normalized = Int[channels...]
    isempty(normalized) && error("At least one channel must be provided")
    any(ch -> ch < 1, normalized) && error("Channels must be 1-based positive indices")
    return normalized
end

function _normalize_pairs(
    pairs::Union{AbstractVector{<:Tuple}, Nothing},
)::Union{Vector{Tuple{Int, Int}}, Nothing}
    pairs === nothing && return nothing
    normalized = Tuple{Int, Int}[]
    for pair in pairs
        length(pair) == 2 || error("Channel pairs must contain exactly two entries")
        first_idx = Int(pair[1])
        second_idx = Int(pair[2])
        first_idx >= 1 || error("Channel pairs must use 1-based positive indices")
        second_idx >= 1 || error("Channel pairs must use 1-based positive indices")
        push!(normalized, (first_idx, second_idx))
    end
    return normalized
end

function _append_passthrough_value!(args::Vector{String}, value)
    if value isa AbstractVector
        for item in value
            push!(args, string(item))
        end
        return args
    end
    push!(args, string(value))
    return args
end

function _append_passthrough!(args::Vector{String}, flag::String, value)
    value === nothing && return args
    push!(args, flag)
    return _append_passthrough_value!(args, value)
end

function _normalize_passthrough_int_list(name::String, value)::Union{Vector{Int}, Nothing}
    value === nothing && return nothing
    if !(value isa AbstractVector || value isa Tuple)
        error("`$name` must be a list of integers")
    end
    normalized = Int[]
    for item in value
        item isa Integer || error("`$name` must contain only integers")
        push!(normalized, Int(item))
    end
    return normalized
end

function _build_passthrough_args(;
    tau_file::Union{AbstractString, Nothing}=nothing,
    tau2=nothing,
    model2=nothing,
    WL_CT::Union{Integer, Nothing}=nothing,
    WS_CT::Union{Integer, Nothing}=nothing,
    no_norm::Bool=false,
    WN_list=nothing,
)::Vector{String}
    args = String[]
    normalized_tau2 = _normalize_passthrough_int_list("tau2", tau2)
    normalized_model2 = _normalize_passthrough_int_list("model2", model2)
    normalized_WN_list = _normalize_passthrough_int_list("WN_list", WN_list)

    _append_passthrough!(args, "-TAU_file", tau_file)
    _append_passthrough!(args, "-TAU2", normalized_tau2)
    _append_passthrough!(args, "-MODEL2", normalized_model2)
    _append_passthrough!(args, "-WL_CT", WL_CT)
    _append_passthrough!(args, "-WS_CT", WS_CT)
    no_norm && push!(args, "-NoNorm")
    _append_passthrough!(args, "-WN_list", normalized_WN_list)
    return args
end

function _normalize_select(
    select::Union{AbstractVector{<:Integer}, Nothing},
)::Union{Vector{Int}, Nothing}
    select === nothing && return nothing
    normalized = Int[select...]
    length(normalized) == SELECT_MASK_SIZE || error(
        "`select` must have $(SELECT_MASK_SIZE) entries",
    )
    all(bit -> bit in (0, 1), normalized) || error("`select` entries must be 0 or 1")
    any(!iszero, normalized) || error("`select` must enable at least one variant")
    return normalized
end

function _resolve_variants(
    variants::Union{AbstractVector{<:AbstractString}, Nothing},
    select::Union{Vector{Int}, Nothing},
)::Vector{String}
    if select !== nothing
        resolved = parse_select_mask(select)
        isempty(resolved) && error("`select` must enable at least one non-reserved variant")
        return resolved
    end

    variants === nothing && error("Provide `variants` or `select`")
    normalized = String[variants...]
    isempty(normalized) && error("At least one variant must be provided")
    return normalized
end

function _resolve_derivative_points(
    derivative_points::Union{Int, Nothing},
    dm::Union{Int, Nothing},
)::Int
    provided = filter(!isnothing, Any[derivative_points, dm])
    if !isempty(provided)
        reference = Int(first(provided))
        for value in provided[2:end]
            Int(value) == reference || error(
                "`derivative_points` and `dm` disagree",
            )
        end
        reference > 0 || error("Derivative points must be positive")
        return reference
    end
    return DDADefaults.DERIVATIVE_POINTS
end

function _validate_custom_model_request(
    model::OptionalModelSpec,
    model_encoding::OptionalModelSpec,
    derivative_points::Union{Int, Nothing},
    dm::Union{Int, Nothing},
    order::Union{Int, Nothing},
)
    explicit_model = model !== nothing || model_encoding !== nothing
    has_derivative_points = !isnothing(derivative_points) || !isnothing(dm)
    if explicit_model && (!has_derivative_points || isnothing(order))
        error(
            "Passing `model` requires explicit `derivative_points`, and `order`",
        )
    end
    return nothing
end

function _resolve_model_terms(
    model::OptionalModelSpec,
    model_encoding::OptionalModelSpec,
    nr_tau::Int,
    order::Int,
)::Vector{Int}
    selected = something(model_encoding, model, copy(DDADefaults.MODEL_PARAMS))
    if selected isa AbstractMatrix
        return model_matrix_to_encoding(
            selected;
            num_delays=nr_tau,
            polynomial_order=order,
        )
    end
    return Int[selected...]
end

function _normalize_sampling_rate(
    sampling_rate::Union{
        Nothing,
        Real,
        Tuple{Real, Real},
        AbstractVector{<:Real},
    },
)::SamplingRate
    sampling_rate === nothing && return nothing

    if sampling_rate isa Real
        return _normalize_sampling_rate_value("Sampling rate", sampling_rate)
    end

    values = sampling_rate isa Tuple ? collect(sampling_rate) : collect(sampling_rate)
    length(values) == 2 || error("Sampling rate must be `nothing`, a scalar, or exactly two numbers")

    first_rate = _normalize_sampling_rate_value("First sampling rate", values[1])
    second_rate = _normalize_sampling_rate_value("Second sampling rate", values[2])
    return (first_rate, second_rate)
end

function _normalize_sampling_rate_value(name::String, value::Real)::Int
    isfinite(Float64(value)) || error("$name must be finite")
    (value isa Integer || isinteger(value)) || error("$name must be integer-valued")
    rate = Int(value)
    rate > 0 || error("$name must be positive")
    return rate
end

function _resolve_tm(
    delays::AbstractVector{<:Integer},
    TM::Union{Int, Nothing},
)::Int
    if TM !== nothing
        TM >= 0 || error("TM must be non-negative")
        return TM
    end
    isempty(delays) && return 0
    return maximum(Int[delays...])
end

function _sampling_rate_scale(
    sampling_rate::SamplingRate,
)::Float64
    sampling_rate === nothing && return 1.0
    sampling_rate isa Integer && return Float64(sampling_rate)
    return Float64(max(sampling_rate[1], sampling_rate[2]))
end

function _sampling_rate_args(sampling_rate::SamplingRate)::Vector{String}
    sampling_rate === nothing && return String[]
    sampling_rate isa Integer && return [string(sampling_rate)]
    return [string(sampling_rate[1]), string(sampling_rate[2])]
end

"""
    DDARequest(file_path, channels, variants; kwargs...)

Create a DDA analysis request.

# Keyword Arguments
- `WL`: Optional analysis window length passed as `-WL`
- `WS`: Optional window step size passed as `-WS`
- `ct_window_length`: CT-specific window length passed as `-WL_CT` when set
- `ct_window_step`: CT-specific window step passed as `-WS_CT` when set
- `delays`: Delay (tau) values, default `$(DEFAULT_DELAYS)`
- `model`: Optional custom model. A vector is passed as `-MODEL` indices; a matrix is converted row-wise from monomial encodings to indices.
- `model_encoding`: Backward-compatible alias for `model`
- `derivative_points::Int=$(DDADefaults.DERIVATIVE_POINTS)`: Value passed to binary `-dm`
- `dm`: Legacy alias for `derivative_points`
- `order::Int=$(DDADefaults.POLYNOMIAL_ORDER)`: Polynomial order. Required when passing a custom `model`
- `nr_tau::Int=$(DDADefaults.NUM_TAU)`: Number of tau values
- `time_range`: Optional `(start, stop)` in samples
- `ct_pairs`: CT channel pairs (1-based)
- `cd_pairs`: CD directed pairs (1-based)
- `select`: Optional explicit `-SELECT` mask. When passed, it overrides `variants`
- `sampling_rate`: Optional `-SR` value. A scalar emits `-SR N`; a tuple emits `-SR N1 N2`. Defaults to `$(DDADefaults.SAMPLING_RATE)`
- `TM`: Optional value used only to compute the derived `t` axis. Defaults to `max(delays)`
- `out_fn`: Optional output base passed to `-OUT_FN`
- `tau_file`, `tau2`, `model2`, `WL_CT`, `WS_CT`, `no_norm`, `WN_list`: Raw binary passthrough arguments
"""
function DDARequest(
    file_path::AbstractString,
    channels::AbstractVector{<:Integer},
    variants::Union{AbstractVector{<:AbstractString}, Nothing}=nothing;
    WL::Union{Int, Nothing}=DDADefaults.WL,
    WS::Union{Int, Nothing}=DDADefaults.WS,
    ct_window_length::Union{Int, Nothing}=nothing,
    ct_window_step::Union{Int, Nothing}=nothing,
    delays::Vector{Int}=collect(DEFAULT_DELAYS),
    model::OptionalModelSpec=nothing,
    model_encoding::OptionalModelSpec=nothing,
    derivative_points::Union{Int, Nothing}=nothing,
    dm::Union{Int, Nothing}=nothing,
    order::Union{Int, Nothing}=nothing,
    nr_tau::Int=DDADefaults.NUM_TAU,
    time_range::Union{Tuple{Real, Real}, Nothing}=nothing,
    ct_pairs::Union{AbstractVector{<:Tuple}, Nothing}=nothing,
    cd_pairs::Union{AbstractVector{<:Tuple}, Nothing}=nothing,
    select::Union{AbstractVector{<:Integer}, Nothing}=nothing,
    sampling_rate::Union{
        Nothing,
        Real,
        Tuple{Real, Real},
        AbstractVector{<:Real},
    }=DDADefaults.SAMPLING_RATE,
    TM::Union{Int, Nothing}=nothing,
    out_fn::Union{AbstractString, Nothing}=nothing,
    tau_file::Union{AbstractString, Nothing}=nothing,
    tau2=nothing,
    model2=nothing,
    WL_CT::Union{Integer, Nothing}=nothing,
    WS_CT::Union{Integer, Nothing}=nothing,
    no_norm::Bool=false,
    WN_list=nothing,
)
    normalized_channels = _normalize_channels(channels)
    normalized_select = _normalize_select(select)
    normalized_variants = _resolve_variants(variants, normalized_select)
    _validate_custom_model_request(
        model,
        model_encoding,
        derivative_points,
        dm,
        order,
    )
    wp = WindowParameters(WL, WS, ct_window_length, ct_window_step)
    dp = DelayParameters(Int[delays...])
    mp = ModelParameters(
        _resolve_derivative_points(derivative_points, dm),
        something(order, DDADefaults.POLYNOMIAL_ORDER),
        nr_tau,
    )
    terms = _resolve_model_terms(model, model_encoding, nr_tau, mp.order)
    tr = time_range === nothing ? nothing : TimeRange(Float64(time_range[1]), Float64(time_range[2]))
    normalized_out_fn = out_fn === nothing ? nothing : expanduser(String(out_fn))
    passthrough_args = _build_passthrough_args(;
        tau_file=tau_file,
        tau2=tau2,
        model2=model2,
        WL_CT=WL_CT,
        WS_CT=WS_CT,
        no_norm=no_norm,
        WN_list=WN_list,
    )
    return DDARequest(
        String(file_path),
        normalized_channels,
        normalized_variants,
        wp,
        dp,
        mp,
        terms,
        tr,
        _normalize_pairs(ct_pairs),
        _normalize_pairs(cd_pairs),
        normalized_select,
        _normalize_sampling_rate(sampling_rate),
        _resolve_tm(dp.delays, TM),
        normalized_out_fn,
        passthrough_args,
    )
end
