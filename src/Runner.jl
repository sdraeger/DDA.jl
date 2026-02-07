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
export DDARunner, run_analysis
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
    window_length::Int
    window_step::Int
    ct_window_length::Union{Int, Nothing}
    ct_window_step::Union{Int, Nothing}
end

WindowParameters(wl::Int, ws::Int) = WindowParameters(wl, ws, nothing, nothing)

"""Delay parameters for DDA analysis."""
struct DelayParameters
    delays::Vector{Int}
end

DelayParameters() = DelayParameters(collect(DEFAULT_DELAYS))

"""Model parameters for DDA."""
struct ModelParameters
    dm::Int
    order::Int
    nr_tau::Int
end

ModelParameters() = ModelParameters(DDADefaults.MODEL_DIMENSION, DDADefaults.POLYNOMIAL_ORDER, DDADefaults.NUM_TAU)

"""
    DDARequest

DDA analysis request parameters.

# Fields
- `file_path`: Path to input file (EDF or ASCII)
- `channels`: Channel indices (0-based)
- `variants`: Variant abbreviations (e.g., ["ST", "SY"])
- `window_params`: Window length and step
- `delay_params`: Delay (tau) values
- `model_params`: Embedding parameters (dm, order, nr_tau)
- `model_encoding`: Model term indices (e.g., [1, 2, 10])
- `time_range`: Optional time range in samples
- `ct_channel_pairs`: Channel pairs for CT (0-based)
- `cd_channel_pairs`: Directed pairs for CD (0-based)
- `sampling_rate`: Optional sampling rate in Hz
"""
struct DDARequest
    file_path::String
    channels::Vector{Int}
    variants::Vector{String}
    window_params::WindowParameters
    delay_params::DelayParameters
    model_params::ModelParameters
    model_encoding::Vector{Int}
    time_range::Union{TimeRange, Nothing}
    ct_channel_pairs::Union{Vector{Tuple{Int, Int}}, Nothing}
    cd_channel_pairs::Union{Vector{Tuple{Int, Int}}, Nothing}
    sampling_rate::Union{Float64, Nothing}
end

"""
    DDARequest(file_path, channels, variants; kwargs...)

Create a DDA analysis request.

# Keyword Arguments
- `window_length::Int=$(DDADefaults.WINDOW_LENGTH)`: Analysis window length
- `window_step::Int=$(DDADefaults.WINDOW_STEP)`: Window step size
- `ct_window_length`: CT-specific window length
- `ct_window_step`: CT-specific window step
- `delays`: Delay (tau) values, default `$(DEFAULT_DELAYS)`
- `model_encoding`: Model term indices, default `$(DDADefaults.MODEL_PARAMS)`
- `dm::Int=$(DDADefaults.MODEL_DIMENSION)`: Embedding dimension
- `order::Int=$(DDADefaults.POLYNOMIAL_ORDER)`: Polynomial order
- `nr_tau::Int=$(DDADefaults.NUM_TAU)`: Number of tau values
- `time_range`: Optional `(start, stop)` in samples
- `ct_pairs`: CT channel pairs (0-based)
- `cd_pairs`: CD directed pairs (0-based)
- `sampling_rate`: Sampling rate in Hz
"""
function DDARequest(
    file_path::String,
    channels::Vector{Int},
    variants::Vector{String};
    window_length::Int=DDADefaults.WINDOW_LENGTH,
    window_step::Int=DDADefaults.WINDOW_STEP,
    ct_window_length::Union{Int, Nothing}=nothing,
    ct_window_step::Union{Int, Nothing}=nothing,
    delays::Vector{Int}=collect(DEFAULT_DELAYS),
    model_encoding::Vector{Int}=copy(DDADefaults.MODEL_PARAMS),
    dm::Int=DDADefaults.MODEL_DIMENSION,
    order::Int=DDADefaults.POLYNOMIAL_ORDER,
    nr_tau::Int=DDADefaults.NUM_TAU,
    time_range::Union{Tuple{Float64, Float64}, Nothing}=nothing,
    ct_pairs::Union{Vector{Tuple{Int, Int}}, Nothing}=nothing,
    cd_pairs::Union{Vector{Tuple{Int, Int}}, Nothing}=nothing,
    sampling_rate::Union{Float64, Nothing}=nothing,
)
    wp = WindowParameters(window_length, window_step, ct_window_length, ct_window_step)
    dp = DelayParameters(delays)
    mp = ModelParameters(dm, order, nr_tau)
    tr = time_range === nothing ? nothing : TimeRange(time_range[1], time_range[2])
    return DDARequest(file_path, channels, variants, wp, dp, mp, model_encoding, tr, ct_pairs, cd_pairs, sampling_rate)
end

# =============================================================================
# STRUCTURED OUTPUT TYPES
# =============================================================================

"""A single timepoint's parsed data for one channel."""
struct StructuredTimepoint
    window_start::Int64
    window_end::Int64
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

"""Result data for a single variant (legacy)."""
struct VariantResultData
    variant_id::String
    variant_name::String
    q_matrix::Matrix{Float64}
    channel_labels::Union{Vector{String}, Nothing}
end

"""Legacy DDA analysis result."""
struct DDAResult
    id::String
    file_path::String
    channels::Vector{String}
    q_matrix::Matrix{Float64}
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
runner = DDARunner()             # auto-discover
runner = DDARunner("/path/to/binary")  # explicit
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
    DDARunner()

Create a DDARunner by auto-discovering the binary.
"""
function DDARunner()
    path = find_binary()
    if path === nothing
        error(
            "DDA binary '$(BINARY_NAME)' not found. " *
            "Set \$$(BINARY_ENV_VAR) or \$$(BINARY_HOME_ENV_VAR), " *
            "or install to one of: $(DEFAULT_BINARY_PATHS)"
        )
    end
    return DDARunner(path)
end

# =============================================================================
# STRUCTURED ANALYSIS (new path — returns per-variant structured data)
# =============================================================================

"""
    run_analysis_structured(runner, request) -> Dict{String, Vector{StructuredChannelData}}

Execute DDA and return fully structured results per variant.
"""
function run_analysis_structured(runner::DDARunner, request::DDARequest)::Dict{String, Vector{StructuredChannelData}}
    if !isfile(request.file_path)
        error("Input file not found: $(request.file_path)")
    end

    analysis_id = string(UUIDs.uuid4())
    output_base = joinpath(tempdir(), "dda_output_$(analysis_id)")
    cmd = build_command(runner, request, output_base)

    try
        run(cmd)
    catch e
        error("DDA execution failed: $e")
    end

    results = Dict{String, Vector{StructuredChannelData}}()
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

    cleanup_temp_files(output_base, request.variants)
    return results
end

"""
    run_analysis_structured(request) -> Dict{String, Vector{StructuredChannelData}}

Execute DDA with auto-discovered binary.
"""
function run_analysis_structured(request::DDARequest)::Dict{String, Vector{StructuredChannelData}}
    runner = DDARunner()
    return run_analysis_structured(runner, request)
end

# =============================================================================
# LEGACY ANALYSIS (backward compat)
# =============================================================================

"""
    run_analysis(runner, request) -> DDAResult

Execute DDA analysis (legacy interface). Prefer `run_st`/`run_ct`/`run_de` instead.
"""
function run_analysis(runner::DDARunner, request::DDARequest)::DDAResult
    if !isfile(request.file_path)
        error("Input file not found: $(request.file_path)")
    end

    analysis_id = string(UUIDs.uuid4())
    output_base = joinpath(tempdir(), "dda_output_$(analysis_id)")
    cmd = build_command(runner, request, output_base)

    try
        run(cmd)
    catch e
        error("DDA execution failed: $e")
    end

    variant_results = parse_results_legacy(request, output_base)
    cleanup_temp_files(output_base, request.variants)

    if isempty(variant_results)
        error("No data extracted from DDA output")
    end

    primary = first(variant_results)
    channel_labels = ["Channel $(ch + 1)" for ch in request.channels]

    return DDAResult(
        analysis_id, request.file_path, channel_labels,
        primary.q_matrix, variant_results,
        request.window_params, request.delay_params,
        string(Dates.now())
    )
end

"""
    run_analysis(request) -> DDAResult

Execute DDA with auto-discovered binary (legacy interface).
"""
function run_analysis(request::DDARequest)::DDAResult
    runner = DDARunner()
    return run_analysis(runner, request)
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

    # Channels (0-based → 1-based)
    push!(args, "-CH_list")
    for ch in request.channels
        push!(args, string(ch + 1))
    end

    # SELECT mask
    mask = generate_select_mask(request.variants)
    push!(args, "-SELECT")
    for bit in mask
        push!(args, string(bit))
    end

    # Model encoding
    push!(args, "-MODEL")
    for p in request.model_encoding
        push!(args, string(p))
    end

    # Delay values
    push!(args, "-TAU")
    for d in request.delay_params.delays
        push!(args, string(d))
    end

    # Window parameters
    wp = request.window_params
    push!(args, "-WL", string(wp.window_length))
    push!(args, "-WS", string(wp.window_step))

    # Model parameters
    mp = request.model_params
    push!(args, "-dm", string(mp.dm))
    push!(args, "-order", string(mp.order))
    push!(args, "-nr_tau", string(mp.nr_tau))

    # CT window parameters if needed
    needs_ct = any(v -> begin
        variant = get_variant_by_abbrev(v)
        variant !== nothing && requires_ct_params(variant)
    end, request.variants)

    if needs_ct || wp.ct_window_length !== nothing
        ct_wl = something(wp.ct_window_length, wp.window_length)
        ct_ws = something(wp.ct_window_step, wp.window_step)
        push!(args, "-WL_CT", string(ct_wl))
        push!(args, "-WS_CT", string(ct_ws))
    end

    # Time bounds
    if request.time_range !== nothing
        tr = request.time_range
        push!(args, "-StartEnd", string(Int(tr.start)), string(Int(tr.stop)))
    end

    # Sampling rate for high-frequency data
    if request.sampling_rate !== nothing && request.sampling_rate > 1000.0
        sr = request.sampling_rate
        push!(args, "-SR", string(Int(sr / 2)), string(Int(sr)))
    end

    if REQUIRES_SHELL_WRAPPER && !Sys.iswindows()
        return `sh $(runner.binary_path) $args`
    else
        return `$(runner.binary_path) $args`
    end
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

    for ch_idx in 0:(num_channels - 1)
        timepoints = StructuredTimepoint[]
        for row in data_rows
            win_start = Int64(row[1])
            win_end = Int64(row[2])
            start_col = 3 + ch_idx * stride  # 1-indexed
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

# =============================================================================
# LEGACY PARSING (kept for backward compat)
# =============================================================================

"""Parse output file extracting only the first coefficient (legacy)."""
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
    q_matrix = Matrix{Float64}(undef, num_channels, num_timepoints)

    for (t, row) in enumerate(data_rows)
        for ch in 1:num_channels
            col_idx = 3 + (ch - 1) * stride
            if col_idx <= length(row)
                q_matrix[ch, t] = row[col_idx]
            end
        end
    end

    return q_matrix
end

function parse_results_legacy(request::DDARequest, output_base::String)::Vector{VariantResultData}
    results = VariantResultData[]
    for variant_abbrev in request.variants
        variant = get_variant_by_abbrev(variant_abbrev)
        variant === nothing && continue

        actual_file = _find_output_file(output_base, variant, variant_abbrev)
        actual_file === nothing && continue

        q_matrix = parse_output_file(actual_file, variant.stride)
        isempty(q_matrix) && continue

        channel_labels = ["Channel $(ch + 1)" for ch in request.channels]
        if size(q_matrix, 1) < length(channel_labels)
            channel_labels = channel_labels[1:size(q_matrix, 1)]
        end

        push!(results, VariantResultData(variant_abbrev, variant.name, q_matrix, channel_labels))
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
