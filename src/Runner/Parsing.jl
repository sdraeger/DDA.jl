# =============================================================================
# STRUCTURED OUTPUT PARSING
# =============================================================================

"""
    parse_output_file_structured(filepath, stride) -> Vector{StructuredChannelData}

Parse a DDA output file extracting ALL coefficients and errors per stride group.

Each row format: `window_start window_end [stride values per channel]*`
For stride=4 (ST/CT): 3 coefficients + 1 error per channel.
For stride=2 (CD): 1 coefficient + 1 error per directed pair.
For stride=1 (DE/SY): 1 value (ergodicity/synchronization measure).
"""
function _read_numeric_rows(filepath::String)::Vector{Vector{Float64}}
    lines = readlines(filepath)
    data_rows = Vector{Vector{Float64}}()

    for line in lines
        stripped = strip(line)
        (isempty(stripped) || startswith(stripped, '#')) && continue
        parts = split(stripped)
        values = tryparse.(Float64, parts)
        all(v -> v !== nothing, values) || continue
        push!(data_rows, Float64[v for v in values])
    end

    return data_rows
end

function parse_output_file_structured(filepath::String, stride::Integer)::Vector{StructuredChannelData}
    data_rows = _read_numeric_rows(filepath)
    isempty(data_rows) && return StructuredChannelData[]

    num_data_cols = length(data_rows[1]) - 2
    if num_data_cols <= 0 || num_data_cols % stride != 0
        @warn "Invalid data format" num_data_cols stride
        return StructuredChannelData[]
    end

    num_channels = div(num_data_cols, stride)
    channels = StructuredChannelData[]

    for ch_idx in 1:num_channels
        timepoints = StructuredTimepoint[]
        for row in data_rows
            win_start = row[1]
            win_end = row[2]
            start_col = 3 + (ch_idx - 1) * stride
            end_col = start_col + stride - 1
            channel_values = row[start_col:end_col]

            if length(channel_values) >= 2
                coeffs = channel_values[1:end-1]
                err = channel_values[end]
            elseif length(channel_values) == 1
                coeffs = Float64[]
                err = channel_values[1]
            else
                coeffs = Float64[]
                err = 0.0
            end
            push!(timepoints, StructuredTimepoint(win_start, win_end, coeffs, err))
        end
        push!(channels, StructuredChannelData(ch_idx, timepoints))
    end

    return channels
end

function _variant_window_spec(
    request::DDARequest,
    variant_abbrev::AbstractString,
)::Union{Tuple{Int, Int}, Nothing}
    wp = request.window_params
    variant = get_variant_by_abbrev(String(variant_abbrev))
    if variant !== nothing && requires_ct_params(variant)
        return nothing
    end
    return wp.WL === nothing || wp.WS === nothing ? nothing : (wp.WL, wp.WS)
end

function _normalized_window_bounds(
    request::DDARequest,
    variant_abbrev::AbstractString,
)::Tuple{Vector{Int64}, Vector{Int64}}
    return _normalized_window_bounds(
        request,
        variant_abbrev,
        0,
    )
end

function _normalized_window_bounds(
    request::DDARequest,
    variant_abbrev::AbstractString,
    n_windows::Integer,
)::Tuple{Vector{Int64}, Vector{Int64}}
    n_windows < 0 && error("n_windows must be non-negative")
    n_windows == 0 && return (Int64[], Int64[])

    spec = _variant_window_spec(request, variant_abbrev)
    spec === nothing && return (Int64[], Int64[])
    (WL, WS) = spec
    first_start = request.time_range === nothing ? 0 : Int(floor(request.time_range.start))
    window_starts = Vector{Int64}(undef, Int(n_windows))
    window_ends = Vector{Int64}(undef, Int(n_windows))

    for window_idx in 1:Int(n_windows)
        window_start = Int64(first_start + (window_idx - 1) * WS)
        window_end = Int64(window_start + WL)
        window_starts[window_idx] = window_start
        window_ends[window_idx] = window_end
    end

    return window_starts, window_ends
end

function _raw_window_bounds(
    channels::Vector{StructuredChannelData},
)::Tuple{Vector{Int64}, Vector{Int64}}
    isempty(channels) && return (Int64[], Int64[])
    starts = Int64[round(Int64, tp.window_start) for tp in channels[1].timepoints]
    ends = Int64[round(Int64, tp.window_end) for tp in channels[1].timepoints]
    return (starts, ends)
end

function _result_window_bounds(
    request::DDARequest,
    variant_abbrev::String,
    channels::Vector{StructuredChannelData},
)::Tuple{Vector{Int64}, Vector{Int64}}
    spec = _variant_window_spec(request, variant_abbrev)
    spec === nothing && return _raw_window_bounds(channels)
    return _normalized_window_bounds(request, variant_abbrev, length(channels[1].timepoints))
end

function _extract_raw_T(
    channels::Vector{StructuredChannelData},
)::Vector{Float64}
    isempty(channels) && return Float64[]
    return [tp.window_start for tp in channels[1].timepoints]
end

function _extract_raw_T_bounds(
    channels::Vector{StructuredChannelData},
)::Matrix{Int64}
    isempty(channels) && return Matrix{Int64}(undef, 0, 2)
    n_windows = length(channels[1].timepoints)
    T = Matrix{Int64}(undef, n_windows, 2)
    for (window_idx, tp) in enumerate(channels[1].timepoints)
        T[window_idx, 1] = _integer_output_index(tp.window_start)
        T[window_idx, 2] = _integer_output_index(tp.window_end)
    end
    return T
end

function _integer_output_index(value::Real)::Int64
    isinteger(value) || error("DDA output T values must be integer-valued, got $value")
    return Int64(value)
end

function _compute_t_axis(
    T::AbstractVector{<:Real},
    derivative_points::Integer,
    tm::Integer,
    sampling_rate::SamplingRate,
)::Vector{Float64}
    denominator = _sampling_rate_scale(sampling_rate)
    return [
        (Float64(raw_T) + 1 + Int(derivative_points) + Int(tm)) / denominator for raw_T in T
    ]
end

function _binary_payload_matrix(
    coefficients::Array{Float64,3},
    errors::Matrix{Float64},
)::Matrix{Float64}
    n_entities, n_windows, n_coeffs = size(coefficients)
    size(errors) == (n_entities, n_windows) || error("error matrix shape must match coefficient entities and windows")

    A = Matrix{Float64}(undef, n_windows, n_entities * (n_coeffs + 1))
    n_windows == 0 && return A

    col = 1
    for entity_idx in 1:n_entities
        for coeff_idx in 1:n_coeffs
            for window_idx in 1:n_windows
                A[window_idx, col] = coefficients[entity_idx, window_idx, coeff_idx]
            end
            col += 1
        end
        for window_idx in 1:n_windows
            A[window_idx, col] = errors[entity_idx, window_idx]
        end
        col += 1
    end

    return A
end

function _coefficient_arrays(
    channels::Vector{StructuredChannelData},
)::Tuple{Array{Float64, 3}, Matrix{Float64}}
    n_entities = length(channels)
    n_windows = isempty(channels) ? 0 : length(channels[1].timepoints)
    n_coefficients = n_windows == 0 ? 0 : length(channels[1].timepoints[1].coefficients)
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

# =============================================================================
# LEGACY PARSING (kept for backward compat)
# =============================================================================

"""Parse output file extracting only the first coefficient (legacy helper)."""
function parse_output_file(filepath::String, stride::Integer)::Matrix{Float64}
    data_rows = _read_numeric_rows(filepath)
    isempty(data_rows) && return Matrix{Float64}(undef, 0, 0)

    num_timepoints = length(data_rows)
    num_data_cols = length(data_rows[1]) - 2

    if num_data_cols <= 0 || num_data_cols % stride != 0
        return Matrix{Float64}(undef, 0, 0)
    end

    num_channels = div(num_data_cols, stride)
    first_coefficients = Matrix{Float64}(undef, num_channels, num_timepoints)

    for (t, row) in enumerate(data_rows)
        for ch in 1:num_channels
            col_idx = 3 + (ch - 1) * stride
            if col_idx <= length(row)
                first_coefficients[ch, t] = row[col_idx]
            end
        end
    end

    return first_coefficients
end

function _channel_labels_for_variant(
    variant::VariantMetadata,
    base_labels::Vector{String},
)::Vector{String}
    if variant.channel_format == Individual
        return copy(base_labels)
    elseif variant.channel_format == Pairs
        labels = String[]
        for i in 1:length(base_labels), j in (i + 1):length(base_labels)
            push!(labels, "$(base_labels[i])-$(base_labels[j])")
        end
        return labels
    elseif variant.channel_format == DirectedPairs
        labels = String[]
        for i in 1:length(base_labels), j in 1:length(base_labels)
            i == j && continue
            push!(labels, "$(base_labels[i])->$(base_labels[j])")
        end
        return labels
    end
    return copy(base_labels)
end

function _pack_variant_result(
    variant_abbrev::String,
    variant::VariantMetadata,
    channels::Vector{StructuredChannelData},
    request::DDARequest,
    channel_labels::Union{Vector{String}, Nothing},
)::VariantResultData
    n_entities = length(channels)
    coefficients, errors = _coefficient_arrays(channels)

    raw_T = _extract_raw_T(channels)
    T = _extract_raw_T_bounds(channels)
    t = _compute_t_axis(
        raw_T,
        request.model_params.derivative_points,
        request.tm,
        request.sampling_rate,
    )
    window_starts, window_ends = _result_window_bounds(request, variant_abbrev, channels)
    A = _binary_payload_matrix(coefficients, errors)

    resolved_labels = channel_labels === nothing ? nothing : begin
        labels = copy(channel_labels)
        if length(labels) > n_entities
            labels = labels[1:n_entities]
        end
        labels
    end

    return VariantResultData(
        variant_abbrev,
        variant.name,
        A,
        coefficients,
        errors,
        T,
        t,
        window_starts,
        window_ends,
        resolved_labels,
    )
end

function parse_results_legacy(request::DDARequest, output_base::String)::Vector{VariantResultData}
    results = VariantResultData[]
    base_labels = request.channels === nothing ? nothing :
                  _resolve_requested_channel_labels(
                      request.file_path,
                      request.channels;
                      fallback_prefix="Channel ",
                  )

    for variant_abbrev in request.variants
        variant = get_variant_by_abbrev(variant_abbrev)
        variant === nothing && continue

        actual_file = _find_output_file(output_base, variant, variant_abbrev)
        actual_file === nothing && continue

        channels = parse_output_file_structured(actual_file, variant.stride)
        isempty(channels) && continue

        variant_labels = base_labels === nothing ?
                         _generic_channel_labels(length(channels)) :
                         _channel_labels_for_variant(variant, base_labels)
        push!(
            results,
            _pack_variant_result(
                variant_abbrev,
                variant,
                channels,
                request,
                variant_labels,
            ),
        )
    end
    return results
end

function _generic_channel_labels(count::Integer)::Vector{String}
    return ["Channel $idx" for idx in 1:Int(count)]
end

# =============================================================================
# HELPERS
# =============================================================================

function _find_output_file(output_base::String, variant::VariantMetadata, abbrev::String)::Union{String, Nothing}
    output_file = "$(output_base)$(variant.output_suffix)"
    isfile(output_file) && return output_file
    legacy_file = "$(output_base)_$(abbrev)"
    isfile(legacy_file) && return legacy_file
    return nothing
end

"""Clean up temporary output files."""
function cleanup_temp_files(output_base::String, variants::Vector{String})
    for variant_abbrev in variants
        variant = get_variant_by_abbrev(variant_abbrev)
        if variant !== nothing
            output_file = "$(output_base)$(variant.output_suffix)"
            isfile(output_file) && rm(output_file, force=true)
        end
        legacy_file = "$(output_base)_$(variant_abbrev)"
        isfile(legacy_file) && rm(legacy_file, force=true)
    end
end
