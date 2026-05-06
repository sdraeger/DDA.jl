#=
DDA Binary Runner Module

Provides functions to execute DDA analysis by running the
run_DDA_AsciiEdf binary and parsing results.
=#

module Runner

import UUIDs
import Dates
using ..Variants
using ..DDADefaults

export DDARequest, DDAResult, VariantResultData
export DDARunner, run_DDA
export StructuredTimepoint, StructuredChannelData
export run_analysis_structured, parse_output_file_structured

# =============================================================================
# TYPES
# =============================================================================

"""Time range for analysis (in samples)."""
struct TimeRange
    start::Float64
    stop::Float64
end

"""Window parameters for DDA analysis."""
struct WindowParameters
    WL::Union{Int, Nothing}
    WS::Union{Int, Nothing}
    ct_window_length::Union{Int, Nothing}
    ct_window_step::Union{Int, Nothing}
end

WindowParameters(WL::Union{Int, Nothing}, WS::Union{Int, Nothing}) = WindowParameters(WL, WS, nothing, nothing)

"""Delay parameters for DDA analysis."""
struct DelayParameters
    delays::Vector{Int}
end

DelayParameters() = DelayParameters(collect(DEFAULT_DELAYS))

"""Model parameters for DDA."""
struct ModelParameters
    derivative_points::Int
    order::Int
    nr_tau::Int
end

ModelParameters() = ModelParameters(
    DDADefaults.DERIVATIVE_POINTS,
    DDADefaults.POLYNOMIAL_ORDER,
    DDADefaults.NUM_TAU,
)

"""
    DDARequest

DDA analysis request parameters.

# Fields
- `file_path`: Path to input file (EDF or ASCII)
- `channels`: Channel indices (1-based, Julia-style)
- `variants`: Variant abbreviations (e.g., ["ST", "SY"])
- `window_params`: Window length and step
- `delay_params`: Delay (tau) values
- `model_params`: Model parameters (`derivative_points`, `order`, `nr_tau`)
- `model_terms`: Model term indices passed to `-MODEL`
- `time_range`: Optional time range in samples
- `ct_channel_pairs`: Channel pairs for CT (1-based)
- `cd_channel_pairs`: Directed pairs for CD (1-based)
- `select`: Optional explicit SELECT mask overriding `variants`
- `sampling_rate`: Optional `-SR` pair `(low, high)`
- `tm`: Optional `TM` value used only to compute the derived `t` axis
- `out_fn`: Optional output base passed to `-OUT_FN`
"""
struct DDARequest
    file_path::String
    channels::Vector{Int}
    variants::Vector{String}
    window_params::WindowParameters
    delay_params::DelayParameters
    model_params::ModelParameters
    model_terms::Vector{Int}
    time_range::Union{TimeRange, Nothing}
    ct_channel_pairs::Union{Vector{Tuple{Int, Int}}, Nothing}
    cd_channel_pairs::Union{Vector{Tuple{Int, Int}}, Nothing}
    select::Union{Vector{Int}, Nothing}
    sampling_rate::Union{Tuple{Int, Int}, Nothing}
    tm::Int
    out_fn::Union{String, Nothing}
    passthrough_args::Vector{String}
end

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
    WL_ct::Union{Integer, Nothing}=nothing,
    WS_ct::Union{Integer, Nothing}=nothing,
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
    _append_passthrough!(args, "-WL_CT", WL_ct)
    _append_passthrough!(args, "-WS_CT", WS_ct)
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
    model::Union{AbstractVector{<:Integer}, Nothing},
    model_encoding::Union{AbstractVector{<:Integer}, Nothing},
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

function _normalize_sampling_rate(
    sampling_rate::Union{
        Nothing,
        Real,
        Tuple{Real, Real},
        AbstractVector{<:Real},
    },
)::Union{Tuple{Int, Int}, Nothing}
    sampling_rate === nothing && return nothing

    if sampling_rate isa Real
        upper = Int(sampling_rate)
        upper > 0 || error("Sampling rate must be positive")
        return (Int(floor(upper / 2)), upper)
    end

    values = sampling_rate isa Tuple ? collect(sampling_rate) : collect(sampling_rate)
    length(values) == 2 || error("Sampling rate must be `nothing`, a scalar, or exactly two numbers")

    low = Int(values[1])
    high = Int(values[2])
    low > 0 || error("Sampling rate lower bound must be positive")
    high > 0 || error("Sampling rate upper bound must be positive")
    low <= high || error("Sampling rate lower bound must be <= upper bound")
    return (low, high)
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
    sampling_rate::Union{Tuple{Int, Int}, Nothing},
)::Float64
    sampling_rate === nothing && return 1.0
    return Float64(max(sampling_rate[1], sampling_rate[2]))
end

function _should_pass_sampling_rate(
    sampling_rate::Union{Tuple{Int, Int}, Nothing},
)::Bool
    sampling_rate === nothing && return false
    return sampling_rate[1] != sampling_rate[2]
end

function _fallback_channel_label(channel::Integer, fallback_prefix::String)::String
    return string(fallback_prefix, Int(channel))
end

function _sanitize_channel_label(label::AbstractString)::String
    cleaned = replace(String(label), '\0' => ' ')
    cleaned = strip(cleaned)
    isempty(cleaned) && return ""
    return strip(cleaned, ['"', '\''])
end

function _read_edf_channel_labels(file_path::AbstractString)::Union{Vector{String}, Nothing}
    open(file_path, "r") do io
        fixed_header = read(io, 256)
        length(fixed_header) == 256 || return nothing

        signal_count = tryparse(Int, strip(String(fixed_header[253:256])))
        signal_count === nothing && return nothing
        signal_count > 0 || return nothing

        labels = String[]
        for _ in 1:signal_count
            field = read(io, 16)
            length(field) == 16 || return nothing
            push!(labels, _sanitize_channel_label(String(field)))
        end
        return labels
    end
end

function _split_ascii_fields(line::AbstractString)::Vector{String}
    stripped = strip(replace(line, '\ufeff' => ' '))
    isempty(stripped) && return String[]

    if occursin('\t', stripped)
        parts = split(stripped, '\t'; keepempty=true)
    elseif occursin(',', stripped)
        parts = split(stripped, ','; keepempty=true)
    else
        parts = split(stripped)
    end

    return [_sanitize_channel_label(part) for part in parts]
end

function _is_numeric_field(field::AbstractString)::Bool
    stripped = strip(field)
    isempty(stripped) && return false
    return tryparse(Float64, stripped) !== nothing
end

function _read_ascii_channel_labels(file_path::AbstractString)::Union{Vector{String}, Nothing}
    return open(file_path, "r") do io
        for line in eachline(io)
            stripped = strip(replace(line, '\ufeff' => ' '))
            (isempty(stripped) || startswith(stripped, '#')) && continue

            fields = _split_ascii_fields(stripped)
            isempty(fields) && continue

            all(_is_numeric_field, fields) && return nothing
            return fields
        end

        return nothing
    end
end

function _infer_input_channel_labels(file_path::AbstractString)::Union{Vector{String}, Nothing}
    isfile(file_path) || return nothing

    try
        ext = lowercase(splitext(String(file_path))[2])
        labels = ext == ".edf" ? _read_edf_channel_labels(file_path) : _read_ascii_channel_labels(file_path)
        labels === nothing && return nothing
        any(!isempty, labels) || return nothing
        return labels
    catch
        return nothing
    end
end

function _resolve_requested_channel_labels(
    file_path::AbstractString,
    channels::AbstractVector{<:Integer};
    fallback_prefix::String="Channel ",
)::Vector{String}
    inferred = _infer_input_channel_labels(file_path)
    resolved = String[]

    for channel in channels
        idx = Int(channel)
        if inferred !== nothing && idx <= length(inferred)
            label = _sanitize_channel_label(inferred[idx])
            if !isempty(label)
                push!(resolved, label)
                continue
            end
        end
        push!(resolved, _fallback_channel_label(idx, fallback_prefix))
    end

    return resolved
end

"""
    DDARequest(file_path, channels, variants; kwargs...)

Create a DDA analysis request.

# Keyword Arguments
- `WL`: Optional analysis window length passed as `-WL`
- `WS`: Optional window step size passed as `-WS`
- `ct_window_length`: CT-specific window length
- `ct_window_step`: CT-specific window step
- `delays`: Delay (tau) values, default `$(DEFAULT_DELAYS)`
- `model`: Optional custom model term indices passed to `-MODEL`
- `model_encoding`: Backward-compatible alias for `model`
- `derivative_points::Int=$(DDADefaults.DERIVATIVE_POINTS)`: Value passed to binary `-dm`
- `dm`: Legacy alias for `derivative_points`
- `order::Int=$(DDADefaults.POLYNOMIAL_ORDER)`: Polynomial order. Required when passing a custom `model`
- `nr_tau::Int=$(DDADefaults.NUM_TAU)`: Number of tau values
- `time_range`: Optional `(start, stop)` in samples
- `ct_pairs`: CT channel pairs (1-based)
- `cd_pairs`: CD directed pairs (1-based)
- `select`: Optional explicit `-SELECT` mask. When passed, it overrides `variants`
- `sampling_rate`: Optional `-SR` pair. Defaults to `$(DDADefaults.SAMPLING_RATE)`
- `TM`: Optional value used only to compute the derived `t` axis. Defaults to `max(delays)`
- `out_fn`: Optional output base passed to `-OUT_FN`
- `tau_file`, `tau2`, `model2`, `WL_ct`, `WS_ct`, `no_norm`, `WN_list`: Raw binary passthrough arguments
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
    model::Union{Vector{Int}, Nothing}=nothing,
    model_encoding::Union{Vector{Int}, Nothing}=nothing,
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
    WL_ct::Union{Integer, Nothing}=nothing,
    WS_ct::Union{Integer, Nothing}=nothing,
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
    terms = Int[something(model_encoding, model, copy(DDADefaults.MODEL_PARAMS))...]
    tr = time_range === nothing ? nothing : TimeRange(Float64(time_range[1]), Float64(time_range[2]))
    normalized_out_fn = out_fn === nothing ? nothing : expanduser(String(out_fn))
    passthrough_args = _build_passthrough_args(;
        tau_file=tau_file,
        tau2=tau2,
        model2=model2,
        WL_ct=WL_ct,
        WS_ct=WS_ct,
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

# =============================================================================
# STRUCTURED OUTPUT TYPES
# =============================================================================

"""A single timepoint's parsed data for one channel."""
struct StructuredTimepoint
    window_start::Float64
    window_end::Float64
    coefficients::Vector{Float64}
    error::Float64
end

"""All timepoints for a single channel/pair."""
struct StructuredChannelData
    channel_index::Int
    timepoints::Vector{StructuredTimepoint}
end

# =============================================================================
# LEGACY RESULT TYPES (backward compat)
# =============================================================================

"""
Result data for a single variant.

`T` contains the first two integer columns emitted by the binary. `A` contains every
remaining column emitted by the binary, preserving output-file order.
"""
struct VariantResultData
    variant_id::String
    variant_name::String
    A::Matrix{Float64}
    coefficients::Array{Float64,3}
    errors::Matrix{Float64}
    T::Matrix{Int64}
    t::Vector{Float64}
    window_starts::Vector{Int64}
    window_ends::Vector{Int64}
    channel_labels::Union{Vector{String}, Nothing}
end

"""DDA analysis result. `T` and `A` mirror the primary variant's raw partition."""
struct DDAResult
    id::String
    file_path::String
    channels::Vector{String}
    T::Matrix{Int64}
    A::Matrix{Float64}
    variant_results::Vector{VariantResultData}
    window_params::WindowParameters
    delay_params::DelayParameters
    created_at::String
end

# =============================================================================
# RUNNER
# =============================================================================

"""
    DDARunner

Handles execution of the run_DDA_AsciiEdf binary.

# Example
```julia
runner = DDARunner()  # auto-discover
runner = DDARunner("/path/to/run_DDA_AsciiEdf")
```
"""
struct DDARunner
    binary_path::String

    function DDARunner(binary_path::AbstractString)
        path = expanduser(binary_path)
        if !isfile(path)
            error("DDA binary not found: $path")
        end
        new(path)
    end
end

"""
    DDARunner(; binary_path=nothing)

Create a runner by resolving the DDA binary from an explicit binary path.
Calling `DDARunner()` with no arguments uses auto-discovery.
"""
function DDARunner(;
    binary_path::Union{AbstractString, Nothing}=nothing,
)
    return DDARunner(require_binary(binary_path))
end

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

"""Internal structured execution helper used by the keyword-only public API."""
function _run_analysis_structured(runner::DDARunner, request::DDARequest)::Dict{String, Vector{StructuredChannelData}}
    if !isfile(request.file_path)
        error("Input file not found: $(request.file_path)")
    end

    output_base, cleanup_output = _resolve_output_base(request)
    cmd = build_command(runner, request, output_base)

    results = Dict{String, Vector{StructuredChannelData}}()
    try
        try
            run(cmd)
        catch e
            error("DDA execution failed: $e")
        end

        for variant_abbrev in request.variants
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
    finally
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
    cmd = build_command(runner, request, output_base)
    variant_results = VariantResultData[]

    try
        try
            run(cmd)
        catch e
            error("DDA execution failed: $e")
        end
        variant_results = parse_results_legacy(request, output_base)
    finally
        cleanup_output && cleanup_temp_files(output_base, request.variants)
    end

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
        primary.T, primary.A, variant_results,
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

# =============================================================================
# COMMAND BUILDING
# =============================================================================

"""Build the DDA CLI command."""
function build_command(runner::DDARunner, request::DDARequest, output_base::String)::Cmd
    args = String[]

    # File type flag
    ext = lowercase(splitext(request.file_path)[2])
    push!(args, ext == ".edf" ? "-EDF" : "-ASCII")

    # Input/output
    push!(args, "-DATA_FN", request.file_path)
    push!(args, "-OUT_FN", output_base)

    # Channels are already 1-based for Julia callers and for the binary.
    push!(args, "-CH_list")
    for ch in request.channels
        push!(args, string(ch))
    end

    # SELECT mask
    mask = request.select === nothing ? generate_select_mask(request.variants) : request.select
    push!(args, "-SELECT")
    for bit in mask
        push!(args, string(bit))
    end

    # Model encoding
    push!(args, "-MODEL")
    for p in request.model_terms
        push!(args, string(p))
    end

    # Delay values
    push!(args, "-TAU")
    for d in request.delay_params.delays
        push!(args, string(d))
    end

    # Window parameters. `-WLms` and `-WSms` are special binary flags and are
    # intentionally not emitted by this wrapper.
    wp = request.window_params
    if wp.WL !== nothing
        push!(args, "-WL", string(wp.WL))
    end
    if wp.WS !== nothing
        push!(args, "-WS", string(wp.WS))
    end

    # Model parameters
    mp = request.model_params
    push!(args, "-dm", string(mp.derivative_points))
    push!(args, "-order", string(mp.order))
    push!(args, "-nr_tau", string(mp.nr_tau))

    # CT window parameters if needed
    needs_ct = any(v -> begin
        variant = get_variant_by_abbrev(v)
        variant !== nothing && requires_ct_params(variant)
    end, request.variants)

    if needs_ct || wp.ct_window_length !== nothing || wp.ct_window_step !== nothing
        ct_wl = wp.ct_window_length === nothing ? wp.WL : wp.ct_window_length
        ct_ws = wp.ct_window_step === nothing ? wp.WS : wp.ct_window_step
        if ct_wl !== nothing && ct_ws !== nothing
            push!(args, "-WL_CT", string(ct_wl))
            push!(args, "-WS_CT", string(ct_ws))
        end
    end

    # Time bounds
    if request.time_range !== nothing
        tr = request.time_range
        push!(args, "-StartEnd", string(Int(tr.start)), string(Int(tr.stop)))
    end

    # Sampling rate pair passed directly to -SR unless it is metadata-only `(N, N)`.
    if _should_pass_sampling_rate(request.sampling_rate)
        push!(args, "-SR", string(request.sampling_rate[1]), string(request.sampling_rate[2]))
    end

    append!(args, request.passthrough_args)

    if REQUIRES_SHELL_WRAPPER && !Sys.iswindows()
        return `sh $(runner.binary_path) $args`
    else
        return `$(runner.binary_path) $args`
    end
end

"""Render the generated shell call as a plain string for diagnostics and tests."""
function build_command_string(
    runner::DDARunner,
    request::DDARequest,
    output_base::String,
)::String
    return join(collect(build_command(runner, request, output_base)), " ")
end

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
function parse_output_file_structured(filepath::String, stride::Integer)::Vector{StructuredChannelData}
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
        WL = wp.ct_window_length === nothing ? wp.WL : wp.ct_window_length
        WS = wp.ct_window_step === nothing ? wp.WS : wp.ct_window_step
        return WL === nothing || WS === nothing ? nothing : (WL, WS)
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
    sampling_rate::Union{Tuple{Int, Int}, Nothing},
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

# =============================================================================
# LEGACY PARSING (kept for backward compat)
# =============================================================================

"""Parse output file extracting only the first coefficient (legacy helper)."""
function parse_output_file(filepath::String, stride::Integer)::Matrix{Float64}
    lines = readlines(filepath)
    isempty(lines) && return Matrix{Float64}(undef, 0, 0)

    data_rows = Vector{Vector{Float64}}()
    for line in lines
        stripped = strip(line)
        (isempty(stripped) || startswith(stripped, '#')) && continue
        parts = split(stripped)
        values = tryparse.(Float64, parts)
        all(v -> v !== nothing, values) || continue
        push!(data_rows, Float64[v for v in values])
    end

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
    n_windows = isempty(channels) ? 0 : length(channels[1].timepoints)
    n_coeffs = n_windows == 0 ? 0 : length(channels[1].timepoints[1].coefficients)

    coefficients = Array{Float64,3}(undef, n_entities, n_windows, n_coeffs)
    errors = Matrix{Float64}(undef, n_entities, n_windows)

    for (entity_idx, channel_data) in enumerate(channels)
        for (window_idx, tp) in enumerate(channel_data.timepoints)
            for (coeff_idx, coeff) in enumerate(tp.coefficients)
                coefficients[entity_idx, window_idx, coeff_idx] = coeff
            end
            errors[entity_idx, window_idx] = tp.error
        end
    end

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
    base_labels = _resolve_requested_channel_labels(
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

        variant_labels = _channel_labels_for_variant(variant, base_labels)
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

end # module Runner
