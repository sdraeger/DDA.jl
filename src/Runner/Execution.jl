function _resolve_output_base(request::DDARequest)::Tuple{String, Bool}
    if request.out_fn !== nothing
        mkpath(dirname(request.out_fn))
        return (request.out_fn, false)
    end
    analysis_id = string(UUIDs.uuid4())
    return (joinpath(tempdir(), "dda_output_$(analysis_id)"), true)
end

# =============================================================================
# STRUCTURED ANALYSIS (new path — returns per-variant structured data)
# =============================================================================

function _request_for_variants(
    request::DDARequest,
    channels::Vector{Int},
    variants::Vector{String},
)::DDARequest
    return DDARequest(
        request.file_path,
        channels,
        variants,
        request.window_params,
        request.delay_params,
        request.model_params,
        request.model_terms,
        request.time_range,
        request.ct_channel_pairs,
        request.cd_channel_pairs,
        nothing,
        request.sampling_rate,
        request.tm,
        request.out_fn,
        request.passthrough_args,
    )
end

function _run_command(runner::DDARunner, request::DDARequest, output_base::String)
    cmd = build_command(runner, request, output_base)
    try
        run(cmd)
    catch e
        error("DDA execution failed: $e")
    end
    return nothing
end

function _logical_command_parts(
    runner::DDARunner,
    request::DDARequest,
    output_base::String,
)::Vector{String}
    parts = collect(build_command(runner, request, output_base))
    if REQUIRES_SHELL_WRAPPER && !Sys.iswindows() && length(parts) >= 2 && parts[1] == SHELL_COMMAND
        return parts[2:end]
    end
    return parts
end

function _write_logical_command_info(
    runner::DDARunner,
    request::DDARequest,
    output_base::String,
)
    open("$(output_base).info", "w") do io
        println(io, join(_logical_command_parts(runner, request, output_base), " "))
    end
    return nothing
end

function _info_file_missing_or_empty(output_base::String)::Bool
    info_file = "$(output_base).info"
    !isfile(info_file) && return true
    return isempty(strip(read(info_file, String)))
end

function _ensure_logical_command_info(
    runner::DDARunner,
    request::DDARequest,
    output_base::String;
    overwrite::Bool=false,
)
    if overwrite || _info_file_missing_or_empty(output_base)
        _write_logical_command_info(runner, request, output_base)
    end
    return nothing
end

function _finalize_output_info(
    runner::DDARunner,
    request::DDARequest,
    output_base::String,
    cleanup_output::Bool,
)
    cleanup_output && return nothing
    return _ensure_logical_command_info(
        runner,
        request,
        output_base;
        overwrite=_has_ct_variant(request),
    )
end

function _parse_structured_outputs(
    request::DDARequest,
    output_base::String,
    variants::Vector{String},
)::Dict{String, Vector{StructuredChannelData}}
    results = Dict{String, Vector{StructuredChannelData}}()
    for variant_abbrev in variants
        variant = get_variant_by_abbrev(variant_abbrev)
        variant === nothing && continue

        actual_file = _find_output_file(output_base, variant, variant_abbrev)
        actual_file === nothing && continue

        channels = parse_output_file_structured(actual_file, variant.stride)
        if !isempty(channels)
            results[variant_abbrev] = channels
        end
    end
    return results
end

function _non_ct_variants(request::DDARequest)::Vector{String}
    return String[v for v in request.variants if v != "CT"]
end

function _has_ct_variant(request::DDARequest)::Bool
    return any(==("CT"), request.variants)
end

function _ct_pairs(request::DDARequest)::Vector{Tuple{Int, Int}}
    if request.ct_channel_pairs !== nothing && !isempty(request.ct_channel_pairs)
        return copy(request.ct_channel_pairs)
    end

    pairs = Tuple{Int, Int}[]
    for i in 1:length(request.channels), j in (i + 1):length(request.channels)
        push!(pairs, (request.channels[i], request.channels[j]))
    end
    isempty(pairs) && error("CT analysis requires at least two channels or an explicit `ct_pairs` list")
    return pairs
end

function _passthrough_int_value(args::Vector{String}, flag::String)::Union{Int, Nothing}
    idx = findfirst(==(flag), args)
    idx === nothing && return nothing
    idx < length(args) || error("`$flag` requires an integer value")
    value = tryparse(Int, args[idx + 1])
    value === nothing && error("`$flag` requires an integer value")
    return value
end

function _merged_ct_value(
    name::String,
    window_value::Union{Int, Nothing},
    passthrough_value::Union{Int, Nothing},
)::Union{Int, Nothing}
    if window_value !== nothing && passthrough_value !== nothing && window_value != passthrough_value
        error("Conflicting `$name` values: $window_value and $passthrough_value")
    end
    return window_value !== nothing ? window_value : passthrough_value
end

function _validate_pairwise_ct_window_args(request::DDARequest)
    wl_ct = _merged_ct_value(
        "WL_CT",
        request.window_params.ct_window_length,
        _passthrough_int_value(request.passthrough_args, "-WL_CT"),
    )
    ws_ct = _merged_ct_value(
        "WS_CT",
        request.window_params.ct_window_step,
        _passthrough_int_value(request.passthrough_args, "-WS_CT"),
    )

    if (wl_ct === nothing) != (ws_ct === nothing)
        error("Pairwise CT requires both `WL_CT` and `WS_CT`, or neither")
    end
    if wl_ct !== nothing && (wl_ct != 2 || ws_ct != 2)
        error(
            "Pairwise CT uses channel-group parameters `WL_CT=2, WS_CT=2`; " *
            "omit them or pass exactly 2/2. Got WL_CT=$wl_ct, WS_CT=$ws_ct.",
        )
    end
    return nothing
end

function _label_map_for_request(request::DDARequest)::Dict{Int, String}
    labels = _resolve_requested_channel_labels(
        request.file_path,
        request.channels;
        fallback_prefix="Channel ",
    )
    return Dict(channel => labels[idx] for (idx, channel) in enumerate(request.channels))
end

function _ct_pair_label(pair::Tuple{Int, Int}, labels::Dict{Int, String})::String
    left = get(labels, pair[1], _fallback_channel_label(pair[1], "Channel "))
    right = get(labels, pair[2], _fallback_channel_label(pair[2], "Channel "))
    return "$left-$right"
end

function _format_output_number(value::Real)::String
    return isinteger(value) ? string(Int(value)) : string(value)
end

function _write_structured_channels(
    filepath::String,
    channels::Vector{StructuredChannelData},
)
    isempty(channels) && return nothing
    n_windows = length(channels[1].timepoints)
    open(filepath, "w") do io
        for window_idx in 1:n_windows
            reference = channels[1].timepoints[window_idx]
            row = String[
                _format_output_number(reference.window_start),
                _format_output_number(reference.window_end),
            ]
            for channel_data in channels
                length(channel_data.timepoints) == n_windows || error("CT pair outputs have inconsistent window counts")
                tp = channel_data.timepoints[window_idx]
                append!(row, string.(tp.coefficients))
                push!(row, string(tp.error))
            end
            println(io, join(row, " "))
        end
    end
    return nothing
end

function _run_ct_pairs(
    runner::DDARunner,
    request::DDARequest,
    output_base::String,
)::Tuple{Vector{StructuredChannelData}, Vector{String}}
    _validate_pairwise_ct_window_args(request)
    pair_data = StructuredChannelData[]
    pair_labels = String[]
    labels = _label_map_for_request(request)

    for (pair_idx, pair) in enumerate(_ct_pairs(request))
        pair_request = _request_for_variants(request, [pair[1], pair[2]], ["CT"])
        pair_output_base = "$(output_base)_ct_pair_$(pair_idx)"
        try
            _run_command(runner, pair_request, pair_output_base)
            actual_file = _find_output_file(pair_output_base, CT, "CT")
            actual_file === nothing && error("No CT output file produced for channel pair $(pair)")
            channels = parse_output_file_structured(actual_file, CT.stride)
            isempty(channels) && error("CT output contained no coefficient data for channel pair $(pair)")
            push!(pair_data, channels[1])
            push!(pair_labels, _ct_pair_label(pair, labels))
        finally
            cleanup_temp_files(pair_output_base, ["CT"])
            info_file = "$(pair_output_base).info"
            isfile(info_file) && rm(info_file, force=true)
        end
    end

    _write_structured_channels("$(output_base)$(CT.output_suffix)", pair_data)
    return pair_data, pair_labels
end

"""Internal structured execution helper used by the keyword-only public API."""
function _run_analysis_structured(runner::DDARunner, request::DDARequest)::Dict{String, Vector{StructuredChannelData}}
    if !isfile(request.file_path)
        error("Input file not found: $(request.file_path)")
    end

    output_base, cleanup_output = _resolve_output_base(request)
    results = Dict{String, Vector{StructuredChannelData}}()
    try
        variants = _non_ct_variants(request)
        if !isempty(variants)
            non_ct_request = _request_for_variants(request, request.channels, variants)
            _run_command(runner, non_ct_request, output_base)
            merge!(results, _parse_structured_outputs(non_ct_request, output_base, variants))
        end

        if _has_ct_variant(request)
            ct_channels, _ = _run_ct_pairs(runner, request, output_base)
            results["CT"] = ct_channels
        end
        return results
    finally
        _finalize_output_info(runner, request, output_base, cleanup_output)
        cleanup_output && cleanup_temp_files(output_base, request.variants)
    end
end

"""
    run_analysis_structured(; file_path, channels, flavors, binary_path=nothing, kwargs...)

Execute DDA without constructing a `DDARequest` explicitly.
"""
function run_analysis_structured(;
    request::Union{DDARequest, Nothing}=nothing,
    runner::Union{DDARunner, Nothing}=nothing,
    file_path::Union{AbstractString, Nothing}=nothing,
    channels::Union{AbstractVector{<:Integer}, Nothing}=nothing,
    flavors::Union{AbstractVector{<:AbstractString}, Nothing}=nothing,
    binary_path::Union{AbstractString, Nothing}=nothing,
    kwargs...,
)::Dict{String, Vector{StructuredChannelData}}
    if request !== nothing
        runner_obj = something(runner, DDARunner(; binary_path=binary_path))
        return _run_analysis_structured(runner_obj, request)
    end

    file_path !== nothing || error("`file_path` keyword is required")
    channels !== nothing || error("`channels` keyword is required")
    runner_obj = DDARunner(; binary_path=binary_path)
    request_obj = DDARequest(file_path, channels, flavors; kwargs...)
    return _run_analysis_structured(runner_obj, request_obj)
end

# =============================================================================
# LEGACY ANALYSIS (backward compat)
# =============================================================================

"""Internal legacy execution helper used by the keyword-only public API."""
function _run_DDA(runner::DDARunner, request::DDARequest)::DDAResult
    if !isfile(request.file_path)
        error("Input file not found: $(request.file_path)")
    end

    analysis_id = string(UUIDs.uuid4())
    output_base, cleanup_output = _resolve_output_base(request)

    variant_results_by_id = Dict{String, VariantResultData}()
    try
        variants = _non_ct_variants(request)
        if !isempty(variants)
            non_ct_request = _request_for_variants(request, request.channels, variants)
            _run_command(runner, non_ct_request, output_base)
            for result in parse_results_legacy(non_ct_request, output_base)
                variant_results_by_id[result.variant_id] = result
            end
        end

        if _has_ct_variant(request)
            ct_channels, ct_labels = _run_ct_pairs(runner, request, output_base)
            variant_results_by_id["CT"] = _pack_variant_result(
                "CT",
                CT,
                ct_channels,
                request,
                ct_labels,
            )
        end
    finally
        _finalize_output_info(runner, request, output_base, cleanup_output)
        cleanup_output && cleanup_temp_files(output_base, request.variants)
    end

    variant_results = VariantResultData[
        variant_results_by_id[variant_id] for variant_id in request.variants
        if haskey(variant_results_by_id, variant_id)
    ]

    if isempty(variant_results)
        error("No data extracted from DDA output")
    end

    primary = first(variant_results)
    channel_labels = _resolve_requested_channel_labels(
        request.file_path,
        request.channels;
        fallback_prefix="Channel ",
    )

    return DDAResult(
        analysis_id, request.file_path, channel_labels,
        primary.T, primary.t, primary.A, variant_results,
        request.window_params, request.delay_params,
        string(Dates.now())
    )
end

"""
    run_DDA(; file_path, channels, flavors, binary_path=nothing, kwargs...)

Execute the DDA binary without constructing a `DDARequest` explicitly.
"""
function run_DDA(;
    request::Union{DDARequest, Nothing}=nothing,
    runner::Union{DDARunner, Nothing}=nothing,
    file_path::Union{AbstractString, Nothing}=nothing,
    channels::Union{AbstractVector{<:Integer}, Nothing}=nothing,
    flavors::Union{AbstractVector{<:AbstractString}, Nothing}=nothing,
    binary_path::Union{AbstractString, Nothing}=nothing,
    kwargs...,
)::DDAResult
    if request !== nothing
        runner_obj = something(runner, DDARunner(; binary_path=binary_path))
        return _run_DDA(runner_obj, request)
    end

    file_path !== nothing || error("`file_path` keyword is required")
    channels !== nothing || error("`channels` keyword is required")
    runner_obj = DDARunner(; binary_path=binary_path)
    request_obj = DDARequest(file_path, channels, flavors; kwargs...)
    return _run_DDA(runner_obj, request_obj)
end
