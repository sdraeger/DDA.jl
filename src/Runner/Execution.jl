function _resolve_output_base(request::DDARequest)::Tuple{String, Bool}
    if request.out_fn !== nothing
        mkpath(dirname(request.out_fn))
        return (request.out_fn, false)
    end
    analysis_id = string(UUIDs.uuid4())
    return (joinpath(tempdir(), "dda_output_$(analysis_id)"), true)
end

function _reject_request_keyword_conflict(; file_path, channels, flavors, kwargs...)
    isempty(kwargs) || error(
        "Pass either `request` or individual analysis keywords, not both; " *
        "unexpected keyword(s): $(join(string.(keys(kwargs)), ", "))",
    )
    (file_path === nothing && channels === nothing && flavors === nothing) || error(
        "Pass either `request` or analysis keywords (`file_path`, `channels`, `flavors`), not both",
    )
    return nothing
end

# =============================================================================
# STRUCTURED ANALYSIS (new path — returns per-variant structured data)
# =============================================================================

function _run_command(cmd::Cmd)
    run(cmd)
    return nothing
end

function _logical_command_parts(cmd::Cmd)::Vector{String}
    parts = collect(cmd)
    if REQUIRES_SHELL_WRAPPER && !Sys.iswindows() && length(parts) >= 2 && parts[1] == SHELL_COMMAND
        return parts[2:end]
    end
    return parts
end

function _logical_command_parts(
    runner::DDARunner,
    request::DDARequest,
    output_base::String,
)::Vector{String}
    return _logical_command_parts(build_command(runner, request, output_base))
end

function _write_logical_command_info(cmd::Cmd, output_base::String)
    open("$(output_base).info", "w") do io
        println(io, join(_logical_command_parts(cmd), " "))
    end
    return nothing
end

function _info_file_missing_or_empty(output_base::String)::Bool
    info_file = "$(output_base).info"
    !isfile(info_file) && return true
    return isempty(strip(read(info_file, String)))
end

function _ensure_logical_command_info(cmd::Cmd, output_base::String; overwrite::Bool=false)
    if overwrite || _info_file_missing_or_empty(output_base)
        _write_logical_command_info(cmd, output_base)
    end
    return nothing
end

function _with_DDA_output(operation::Function, runner::DDARunner, request::DDARequest)
    isfile(request.file_path) || error("Input file not found: $(request.file_path)")
    output_base, cleanup_output = _resolve_output_base(request)
    cmd = build_command(runner, request, output_base)
    try
        _run_command(cmd)
        return operation(output_base)
    finally
        cleanup_output || _ensure_logical_command_info(cmd, output_base)
        cleanup_output && cleanup_temp_files(output_base, request.variants)
    end
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

    missing_variants = [v for v in variants if !haskey(results, v)]
    isempty(missing_variants) || error(
        "No DDA output found for requested variant(s): $(join(missing_variants, ", ")). " *
        "The binary may have failed or produced no output for them.",
    )
    return results
end

"""Run the binary once and parse its outputs into packed per-variant data.

This is the single parse path shared by `run_DDA` (legacy result object) and
the typed high-level API (`run_st`/`run_ct`/`run_de`).
"""
function _run_and_parse(runner::DDARunner, request::DDARequest)::Vector{VariantResultData}
    return _with_DDA_output(runner, request) do output_base
        return parse_results_legacy(request, output_base)
    end
end

"""Internal structured execution helper used by the keyword-only public API."""
function _run_analysis_structured(runner::DDARunner, request::DDARequest)::Dict{String, Vector{StructuredChannelData}}
    return _with_DDA_output(runner, request) do output_base
        return _parse_structured_outputs(request, output_base, request.variants)
    end
end

"""
    run_analysis_structured(; file_path, channels=nothing, flavors, binary_path=nothing, kwargs...)

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
        _reject_request_keyword_conflict(; file_path=file_path, channels=channels, flavors=flavors, kwargs...)
        runner_obj = something(runner, DDARunner(; binary_path=binary_path))
        return _run_analysis_structured(runner_obj, request)
    end

    file_path !== nothing || error("`file_path` keyword is required")
    runner_obj = DDARunner(; binary_path=binary_path)
    request_obj = DDARequest(file_path, channels, flavors; kwargs...)
    return _run_analysis_structured(runner_obj, request_obj)
end

# =============================================================================
# LEGACY ANALYSIS (backward compat)
# =============================================================================

"""Internal legacy execution helper used by the keyword-only public API."""
function _run_DDA(runner::DDARunner, request::DDARequest)::DDAResult
    analysis_id = string(UUIDs.uuid4())
    variant_results = _run_and_parse(runner, request)

    if isempty(variant_results)
        error("No data extracted from DDA output")
    end

    primary = first(variant_results)
    channel_labels = request.channels === nothing ? something(primary.channel_labels, String[]) :
                     _resolve_requested_channel_labels(
                         request.file_path,
                         request.channels;
                         fallback_prefix="Channel ",
                     )

    return DDAResult(
        analysis_id, request.file_path, channel_labels,
        primary.T, primary.t, primary.A, variant_results,
        request.window_params, request.delay_params,
        Dates.now()
    )
end

"""Execute the DDA binary without parsing result files into Julia."""
function _run_DDA_command_only(runner::DDARunner, request::DDARequest)::Nothing
    _with_DDA_output(runner, request) do _
        return nothing
    end
    return nothing
end

"""
    run_DDA(; file_path, channels=nothing, flavors, binary_path=nothing, load_results=true, kwargs...)

Execute the DDA binary without constructing a `DDARequest` explicitly. Pass
`load_results=false` to execute the shell command without parsing output files
into Julia; in that mode the function returns `nothing`.
"""
function run_DDA(;
    request::Union{DDARequest, Nothing}=nothing,
    runner::Union{DDARunner, Nothing}=nothing,
    file_path::Union{AbstractString, Nothing}=nothing,
    channels::Union{AbstractVector{<:Integer}, Nothing}=nothing,
    flavors::Union{AbstractVector{<:AbstractString}, Nothing}=nothing,
    binary_path::Union{AbstractString, Nothing}=nothing,
    load_results::Bool=true,
    kwargs...,
)::Union{DDAResult, Nothing}
    if request !== nothing
        _reject_request_keyword_conflict(; file_path=file_path, channels=channels, flavors=flavors, kwargs...)
        runner_obj = something(runner, DDARunner(; binary_path=binary_path))
        if load_results
            return _run_DDA(runner_obj, request)
        end
        return _run_DDA_command_only(runner_obj, request)
    end

    file_path !== nothing || error("`file_path` keyword is required")
    runner_obj = DDARunner(; binary_path=binary_path)
    request_obj = DDARequest(file_path, channels, flavors; kwargs...)
    if load_results
        return _run_DDA(runner_obj, request_obj)
    end
    return _run_DDA_command_only(runner_obj, request_obj)
end
