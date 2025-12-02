#=
DDA Binary Runner Module

Provides high-level functions to execute DDA analysis by running the
run_DDA_AsciiEdf binary and parsing results.
=#

module Runner

import UUIDs
import Dates
using ..Variants

export DDARequest, DDAResult, VariantResultData
export DDARunner, run_analysis, run_analysis!

# =============================================================================
# TYPES
# =============================================================================

"""
Time range for analysis (in samples or seconds depending on context).
"""
struct TimeRange
    start::Float64
    stop::Float64
end

"""
Window parameters for DDA analysis.
"""
struct WindowParameters
    window_length::Int
    window_step::Int
    ct_window_length::Union{Int, Nothing}
    ct_window_step::Union{Int, Nothing}
end

WindowParameters(wl::Int, ws::Int) = WindowParameters(wl, ws, nothing, nothing)

"""
Delay parameters for DDA analysis.
"""
struct DelayParameters
    delays::Vector{Int}
end

DelayParameters() = DelayParameters(DEFAULT_DELAYS)

"""
Model parameters for DDA (expert mode).
"""
struct ModelParameters
    dm::Int      # Embedding dimension (default: 4)
    order::Int   # Polynomial order (default: 4)
    nr_tau::Int  # Number of tau values (default: 2)
end

ModelParameters() = ModelParameters(4, 4, 2)

"""
DDA Analysis Request Parameters.

# Fields
- `file_path`: Path to input file (EDF or ASCII)
- `channels`: Channel indices (0-based) to analyze
- `variants`: Variant abbreviations to run (e.g., ["ST", "SY"])
- `window_params`: Window length and step parameters
- `delay_params`: Delay (tau) values
- `model_params`: Model parameters (dm, order, nr_tau)
- `time_range`: Optional time range for analysis (in samples)
- `ct_channel_pairs`: Channel pairs for CT analysis (0-based)
- `cd_channel_pairs`: Directed channel pairs for CD analysis (0-based)
- `sampling_rate`: Optional sampling rate in Hz
"""
struct DDARequest
    file_path::String
    channels::Vector{Int}
    variants::Vector{String}
    window_params::WindowParameters
    delay_params::DelayParameters
    model_params::ModelParameters
    time_range::Union{TimeRange, Nothing}
    ct_channel_pairs::Union{Vector{Tuple{Int, Int}}, Nothing}
    cd_channel_pairs::Union{Vector{Tuple{Int, Int}}, Nothing}
    sampling_rate::Union{Float64, Nothing}
end

"""
    DDARequest(file_path, channels, variants; kwargs...)

Create a DDA analysis request.

# Arguments
- `file_path::String`: Path to input file
- `channels::Vector{Int}`: Channel indices (0-based)
- `variants::Vector{String}`: Variants to run (e.g., ["ST", "CT"])

# Keyword Arguments
- `window_length::Int=2048`: Analysis window length
- `window_step::Int=1024`: Window step size
- `ct_window_length::Union{Int,Nothing}=nothing`: CT-specific window length
- `ct_window_step::Union{Int,Nothing}=nothing`: CT-specific window step
- `delays::Vector{Int}=DEFAULT_DELAYS`: Delay (tau) values
- `dm::Int=4`: Embedding dimension
- `order::Int=4`: Polynomial order
- `nr_tau::Int=2`: Number of tau values
- `time_range::Union{Tuple{Float64,Float64},Nothing}=nothing`: Time range (samples)
- `ct_pairs::Union{Vector{Tuple{Int,Int}},Nothing}=nothing`: CT channel pairs
- `cd_pairs::Union{Vector{Tuple{Int,Int}},Nothing}=nothing`: CD directed pairs
- `sampling_rate::Union{Float64,Nothing}=nothing`: Sampling rate in Hz

# Example
```julia
request = DDARequest(
    "data.edf",
    [0, 1, 2],  # channels (0-based)
    ["ST", "SY"];
    window_length=2048,
    window_step=1024,
    delays=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
)
```
"""
function DDARequest(
    file_path::String,
    channels::Vector{Int},
    variants::Vector{String};
    window_length::Int=2048,
    window_step::Int=1024,
    ct_window_length::Union{Int, Nothing}=nothing,
    ct_window_step::Union{Int, Nothing}=nothing,
    delays::Vector{Int}=collect(DEFAULT_DELAYS),
    dm::Int=4,
    order::Int=4,
    nr_tau::Int=2,
    time_range::Union{Tuple{Float64, Float64}, Nothing}=nothing,
    ct_pairs::Union{Vector{Tuple{Int, Int}}, Nothing}=nothing,
    cd_pairs::Union{Vector{Tuple{Int, Int}}, Nothing}=nothing,
    sampling_rate::Union{Float64, Nothing}=nothing,
)
    wp = WindowParameters(window_length, window_step, ct_window_length, ct_window_step)
    dp = DelayParameters(delays)
    mp = ModelParameters(dm, order, nr_tau)
    tr = time_range === nothing ? nothing : TimeRange(time_range[1], time_range[2])

    return DDARequest(file_path, channels, variants, wp, dp, mp, tr, ct_pairs, cd_pairs, sampling_rate)
end

"""
Result data for a single variant.
"""
struct VariantResultData
    variant_id::String
    variant_name::String
    q_matrix::Matrix{Float64}
    channel_labels::Union{Vector{String}, Nothing}
end

"""
DDA Analysis Result.

# Fields
- `id`: Unique analysis ID
- `file_path`: Input file path
- `channels`: Channel labels
- `q_matrix`: Primary Q matrix [channels × timepoints]
- `variant_results`: Results for each variant
- `window_params`: Window parameters used
- `delay_params`: Delay parameters used
- `created_at`: Timestamp of analysis
"""
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
DDA Binary Runner.

Handles execution of the run_DDA_AsciiEdf binary.

# Example
```julia
# Auto-discover binary
runner = DDARunner()

# Or specify explicit path
runner = DDARunner("/path/to/run_DDA_AsciiEdf")

result = run_analysis(runner, request)
```
"""
struct DDARunner
    binary_path::String

    # Inner constructor validates binary exists
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

Create a DDARunner by auto-discovering the binary location.

Uses `find_binary()` to locate the DDA binary via:
1. `\$DDA_BINARY_PATH` environment variable
2. `\$DDA_HOME/bin/` directory
3. Default search paths (`~/.local/bin`, `~/bin`, `/usr/local/bin`, `/opt/dda/bin`)

# Returns
- `DDARunner` instance

# Throws
- `ErrorException`: If binary not found

# Examples
```julia
runner = DDARunner()  # Auto-discover
result = run_analysis(runner, request)
```
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

"""
    run_analysis(runner::DDARunner, request::DDARequest) -> DDAResult

Execute DDA analysis with the given request parameters.

# Arguments
- `runner`: DDARunner instance with binary path
- `request`: DDA analysis configuration

# Returns
- `DDAResult` containing the Q matrices and metadata

# Example
```julia
runner = DDARunner()
request = DDARequest("data.edf", [0, 1, 2], ["ST", "SY"])
result = run_analysis(runner, request)
println("Q matrix size: ", size(result.q_matrix))
```
"""
function run_analysis(runner::DDARunner, request::DDARequest)::DDAResult
    # Validate input file
    if !isfile(request.file_path)
        error("Input file not found: $(request.file_path)")
    end

    # Generate unique analysis ID
    analysis_id = string(UUIDs.uuid4())

    # Create temporary output file
    temp_dir = tempdir()
    output_base = joinpath(temp_dir, "dda_output_$(analysis_id)")

    # Build and execute command
    cmd = build_command(runner, request, output_base)

    @info "Executing DDA command" command=cmd

    try
        run(cmd)
        @info "DDA execution completed"
    catch e
        error("DDA execution failed: $e")
    end

    # Parse results
    variant_results = parse_results(request, output_base)

    # Clean up temporary files
    cleanup_temp_files(output_base, request.variants)

    if isempty(variant_results)
        error("No data extracted from DDA output")
    end

    # Use first variant as primary
    primary = first(variant_results)

    # Build channel labels
    channel_labels = ["Channel $(ch + 1)" for ch in request.channels]

    return DDAResult(
        analysis_id,
        request.file_path,
        channel_labels,
        primary.q_matrix,
        variant_results,
        request.window_params,
        request.delay_params,
        string(Dates.now())
    )
end

"""
    run_analysis(request::DDARequest) -> DDAResult

Execute DDA analysis using auto-discovered binary.
"""
function run_analysis(request::DDARequest)::DDAResult
    runner = DDARunner()
    return run_analysis(runner, request)
end

# =============================================================================
# INTERNAL FUNCTIONS
# =============================================================================

"""Build the DDA command with all arguments."""
function build_command(runner::DDARunner, request::DDARequest, output_base::String)::Cmd
    args = String[]

    # File type flag
    ext = lowercase(splitext(request.file_path)[2])
    file_type_flag = ext == ".edf" ? "-EDF" : "-ASCII"
    push!(args, file_type_flag)

    # Input/output files
    push!(args, "-DATA_FN", request.file_path)
    push!(args, "-OUT_FN", output_base)

    # Channel list (convert 0-based to 1-based)
    push!(args, "-CH_list")
    for ch in request.channels
        push!(args, string(ch + 1))
    end

    # Model parameters
    mp = request.model_params
    push!(args, "-dm", string(mp.dm))
    push!(args, "-order", string(mp.order))
    push!(args, "-nr_tau", string(mp.nr_tau))

    # Window parameters
    wp = request.window_params
    push!(args, "-WL", string(wp.window_length))
    push!(args, "-WS", string(wp.window_step))

    # SELECT mask
    mask = generate_select_mask(request.variants)
    push!(args, "-SELECT")
    for bit in mask
        push!(args, string(bit))
    end

    # Model encoding
    push!(args, "-MODEL", "1", "2", "10")

    # CT window parameters if needed
    needs_ct_params = any(v -> begin
        variant = get_variant_by_abbrev(v)
        variant !== nothing && requires_ct_params(variant)
    end, request.variants)

    if needs_ct_params || wp.ct_window_length !== nothing
        ct_wl = something(wp.ct_window_length, 2)
        ct_ws = something(wp.ct_window_step, 2)
        push!(args, "-WL_CT", string(ct_wl))
        push!(args, "-WS_CT", string(ct_ws))
    end

    # Delay values
    push!(args, "-TAU")
    for d in request.delay_params.delays
        push!(args, string(d))
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

    # APE binary needs to run through sh on Unix
    if REQUIRES_SHELL_WRAPPER && !Sys.iswindows()
        return `sh $(runner.binary_path) $args`
    else
        return `$(runner.binary_path) $args`
    end
end

"""Parse DDA output files for all requested variants."""
function parse_results(request::DDARequest, output_base::String)::Vector{VariantResultData}
    results = VariantResultData[]

    for variant_abbrev in request.variants
        variant = get_variant_by_abbrev(variant_abbrev)
        if variant === nothing
            @warn "Unknown variant: $variant_abbrev"
            continue
        end

        # Construct output file path
        output_file = "$(output_base)$(variant.output_suffix)"

        # Also try legacy naming
        legacy_file = "$(output_base)_$(variant_abbrev)"

        actual_file = if isfile(output_file)
            output_file
        elseif isfile(legacy_file)
            legacy_file
        else
            @warn "Output file not found for variant $variant_abbrev" tried=[output_file, legacy_file]
            continue
        end

        @info "Parsing output for variant $variant_abbrev" file=actual_file

        # Parse output file
        q_matrix = parse_output_file(actual_file, variant.stride)

        if isempty(q_matrix)
            @warn "No data extracted for variant $variant_abbrev"
            continue
        end

        # Generate channel labels
        channel_labels = ["Channel $(ch + 1)" for ch in request.channels]
        if size(q_matrix, 1) < length(channel_labels)
            channel_labels = channel_labels[1:size(q_matrix, 1)]
        end

        push!(results, VariantResultData(
            variant_abbrev,
            variant.name,
            q_matrix,
            channel_labels
        ))
    end

    return results
end

"""
Parse a DDA output file into a Q matrix.

# Arguments
- `filepath`: Path to output file
- `stride`: Column stride for this variant (values per channel)

# Returns
Matrix{Float64} with shape [channels × timepoints]
"""
function parse_output_file(filepath::String, stride::UInt8)::Matrix{Float64}
    lines = readlines(filepath)

    if isempty(lines)
        return Matrix{Float64}(undef, 0, 0)
    end

    # Parse all data rows
    data_rows = Vector{Vector{Float64}}()
    for line in lines
        stripped = strip(line)
        if isempty(stripped) || startswith(stripped, '#')
            continue
        end

        parts = split(stripped)
        values = tryparse.(Float64, parts)
        if all(v -> v !== nothing, values)
            push!(data_rows, Float64[v for v in values])
        end
    end

    if isempty(data_rows)
        return Matrix{Float64}(undef, 0, 0)
    end

    # Calculate dimensions
    # Format: window_start window_end [stride values per channel] * num_channels
    num_timepoints = length(data_rows)
    num_data_cols = length(data_rows[1]) - 2  # Exclude window bounds

    if num_data_cols <= 0 || num_data_cols % stride != 0
        @warn "Invalid data format" num_data_cols stride
        return Matrix{Float64}(undef, 0, 0)
    end

    num_channels = div(num_data_cols, stride)

    # Build matrix [channels × timepoints]
    # Extract first coefficient from each stride group
    q_matrix = Matrix{Float64}(undef, num_channels, num_timepoints)

    for (t, row) in enumerate(data_rows)
        for ch in 1:num_channels
            # Skip 2 window cols, then extract first value of each stride group
            col_idx = 3 + (ch - 1) * stride  # 1-indexed
            if col_idx <= length(row)
                q_matrix[ch, t] = row[col_idx]
            end
        end
    end

    return q_matrix
end

"""Clean up temporary output files."""
function cleanup_temp_files(output_base::String, variants::Vector{String})
    for variant_abbrev in variants
        variant = get_variant_by_abbrev(variant_abbrev)
        if variant !== nothing
            output_file = "$(output_base)$(variant.output_suffix)"
            isfile(output_file) && rm(output_file, force=true)
        end
        # Also try legacy naming
        legacy_file = "$(output_base)_$(variant_abbrev)"
        isfile(legacy_file) && rm(legacy_file, force=true)
    end
end

end # module Runner
