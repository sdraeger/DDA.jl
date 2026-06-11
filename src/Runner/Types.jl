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
- `channels`: Optional channel indices (1-based, Julia-style). `nothing` omits `-CH_list`.
- `variants`: Variant abbreviations (e.g., ["ST", "SY"])
- `window_params`: Window length and step
- `delay_params`: Delay (tau) values
- `model_params`: Model parameters (`derivative_points`, `order`, `nr_tau`)
- `model_terms`: Model term indices passed to `-MODEL`
- `time_range`: Optional time range in samples
- `ct_channel_pairs`: Channel pairs for CT (1-based)
- `cd_channel_pairs`: Directed pairs for CD (1-based)
- `select`: Optional explicit SELECT mask overriding `variants`
- `input_format`: Input format used for the binary file-type flag
- `sampling_rate`: Optional `-SR` value. A scalar emits `-SR N`; a tuple emits `-SR N1 N2`
- `tm`: Optional `TM` value used only to compute the derived `t` axis
- `out_fn`: Optional output base passed to `-OUT_FN`
"""
const SamplingRate = Union{Int, Tuple{Int, Int}, Nothing}
const ModelSpec = Union{AbstractVector{<:Integer}, AbstractMatrix{<:Integer}}
const OptionalModelSpec = Union{ModelSpec, Nothing}

struct DDARequest
    file_path::String
    channels::Union{Vector{Int}, Nothing}
    variants::Vector{String}
    window_params::WindowParameters
    delay_params::DelayParameters
    model_params::ModelParameters
    model_terms::Vector{Int}
    time_range::Union{TimeRange, Nothing}
    ct_channel_pairs::Union{Vector{Tuple{Int, Int}}, Nothing}
    cd_channel_pairs::Union{Vector{Tuple{Int, Int}}, Nothing}
    select::Union{Vector{Int}, Nothing}
    input_format::FileType
    sampling_rate::SamplingRate
    tm::Int
    out_fn::Union{String, Nothing}
    passthrough_args::Vector{String}
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
