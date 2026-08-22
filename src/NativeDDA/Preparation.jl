const NATIVE_WINDOW_LENGTH = 200
const NATIVE_WINDOW_STEP = 100

function _validate_samples(samples::AbstractMatrix{<:Real})::Matrix{Float64}
    isempty(samples) && error("DDA samples must not be empty")
    size(samples, 1) > 0 || error("DDA samples must contain at least one row")
    size(samples, 2) > 0 || error("DDA samples must contain at least one channel")
    return Float64.(samples)
end

function _normalize_channel_labels(labels, channel_count::Int)::Vector{String}
    labels === nothing && return ["Ch $channel" for channel in 1:channel_count]
    length(labels) == channel_count || error(
        "Expected $channel_count channel labels, got $(length(labels))",
    )
    return String[string(label) for label in labels]
end

function _normalize_channels(channels, channel_count::Int)::Vector{Int}
    selected = channels === nothing || isempty(channels) ? collect(1:channel_count) : Int[channels...]
    all(channel -> 1 <= channel <= channel_count, selected) || error(
        "Channel indices must be in 1:$channel_count",
    )
    return selected
end

function _normalize_flavor(flavor::AbstractString)::String
    token = replace(uppercase(strip(flavor)), '-' => '_', ' ' => '_')
    aliases = Dict(
        "SINGLE_TIMESERIES" => "ST",
        "CROSS_TIMESERIES" => "CT",
        "CROSS_DYNAMICAL" => "CD",
        "DYNAMICAL_ERGODICITY" => "DE",
        "DELAY_EMBEDDING" => "DE",
        "SYNCHRONIZATION" => "SY",
        "SYNCHRONY" => "SY",
    )
    return get(aliases, token, token)
end

function _build_model_spec(;
    window_length::Int,
    window_step::Int,
    delays,
    model_terms,
    derivative_points::Int,
    order::Int,
    nr_tau::Int,
)::ModelSpec
    derivative_points > 0 || error("derivative_points must be greater than zero")
    nr_tau > 0 || error("nr_tau must be greater than zero")
    window_length > 0 || error("window_length must be greater than zero")
    window_step > 0 || error("window_step must be greater than zero")
    length(delays) >= nr_tau || error(
        "Received $(length(delays)) delays but nr_tau=$nr_tau requires at least $nr_tau",
    )

    selected_delays = Int[delays[index] for index in 1:nr_tau]
    all(>=(0), selected_delays) || error("DDA delays must be non-negative")
    monomials = generate_monomials(nr_tau, order)
    terms = Vector{Vector{Int}}()
    for term in model_terms
        index = Int(term)
        1 <= index <= length(monomials) || error(
            "Model term $index is out of range for the monomial table",
        )
        push!(terms, Int[selected_delays[entry] for entry in monomials[index] if entry != 0])
    end
    isempty(terms) && error("model_terms must contain at least one term")

    return ModelSpec(
        derivative_points,
        window_length,
        window_step,
        maximum(selected_delays),
        terms,
    )
end

function _analysis_bounds(start::Real, stop, row_count::Int)::Tuple{Int,Int}
    start_index = floor(Int, max(Float64(start), 0.0))
    stop_value = stop === nothing ? Float64(row_count - 1) : Float64(stop)
    stop_index = isfinite(stop_value) ? floor(Int, max(stop_value, 0.0)) : row_count - 1
    clamped_start = min(start_index, row_count - 1)
    clamped_stop = max(min(stop_index, row_count - 1), clamped_start)
    return clamped_start, clamped_stop - clamped_start + 1
end

function _prepare_window(
    data::Matrix{Float64},
    bounds_start::Int,
    model::ModelSpec,
    window_index::Int,
    normalization::String,
    nr_exclude::Int,
    derivative_step::Int,
)::PreparedWindow
    marker = model.window_length + model.max_delay + 2model.derivative_points
    first_row = bounds_start + (window_index - 1) * model.window_step + 1
    last_row = first_row + marker - 1
    if last_row <= size(data, 1)
        raw = copy(data[first_row:last_row, :])
    else
        available = copy(data[first_row:end, :])
        filler = isempty(available) ? NaN : available[end, end]
        padding = fill(filler, marker - size(available, 1), size(data, 2))
        raw = vcat(available, padding)
    end

    _apply_nan_runs!(raw, nr_exclude)
    derivative = _central_derivative(raw, model.derivative_points, derivative_step)
    shifted, trimmed = _normalize_window(
        raw,
        derivative,
        model.derivative_points,
        model.max_delay,
        normalization,
    )
    return PreparedWindow(shifted, trimmed, model.max_delay)
end

function _apply_nan_runs!(data::Matrix{Float64}, minimum_run::Int)
    minimum_run <= 0 && return data
    for channel in axes(data, 2)
        run_start = 1
        while run_start <= size(data, 1)
            run_end = run_start
            while run_end < size(data, 1) && data[run_end + 1, channel] == data[run_start, channel]
                run_end += 1
            end
            run_end - run_start + 1 >= minimum_run && (data[run_start:run_end, channel] .= NaN)
            run_start = run_end + 1
        end
    end
    return data
end

function _central_derivative(
    data::Matrix{Float64},
    derivative_points::Int,
    derivative_step::Int,
)::Matrix{Float64}
    rows, channels = size(data)
    rows > 2derivative_points || error(
        "Need more than 2*derivative_points=$(2derivative_points) rows, got $rows",
    )
    step = max(derivative_step, 1)
    stencil_count = div(derivative_points, step)
    stencil_count > 0 || error(
        "Invalid derivative_step=$step for derivative_points=$derivative_points",
    )

    derivative = fill(NaN, channels, rows - 2derivative_points)
    for channel in 1:channels
        for center in (derivative_points + 1):(rows - derivative_points)
            valid = !isnan(data[center, channel])
            value = 0.0
            for stencil in 1:stencil_count
                offset = stencil * step
                plus = data[center + offset, channel]
                minus = data[center - offset, channel]
                (isnan(plus) || isnan(minus)) && (valid = false)
                valid && (value += (plus - minus) / stencil)
            end
            valid && (derivative[channel, center - derivative_points] = value / stencil_count)
        end
    end
    return derivative
end

function _normalize_window(
    raw::Matrix{Float64},
    derivative::Matrix{Float64},
    derivative_points::Int,
    max_delay::Int,
    mode::String,
)::Tuple{Matrix{Float64},Matrix{Float64}}
    shifted_rows = size(raw, 1) - 2derivative_points
    window_length = shifted_rows - max_delay
    window_length >= 0 || error("Window length became negative after max(delay) trim")
    shifted = fill(NaN, shifted_rows, size(raw, 2))
    trimmed = fill(NaN, size(raw, 2), window_length)

    for channel in axes(raw, 2)
        shifted[:, channel] = raw[(derivative_points + 1):(derivative_points + shifted_rows), channel]
        if mode == "raw"
            trimmed[channel, :] = derivative[channel, (max_delay + 1):(max_delay + window_length)]
            continue
        end

        valid = filter(value -> !isnan(value), shifted[:, channel])
        isempty(valid) && continue
        if mode == "minmax"
            offset = minimum(valid)
            scale = maximum(valid) - offset
        elseif mode == "zscore"
            length(valid) >= 2 || continue
            offset = mean(valid)
            scale = std(valid; corrected=true)
        else
            error("Unknown normalization mode $mode")
        end
        (!isfinite(scale) || iszero(scale)) && continue
        shifted[:, channel] = (shifted[:, channel] .- offset) ./ scale
        trimmed[channel, :] = derivative[channel, (max_delay + 1):(max_delay + window_length)] ./ scale
    end
    return shifted, trimmed
end
