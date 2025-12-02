#=
AUTO-GENERATED from DDA Specification v1.0.0
Generated: 2025-12-02

DO NOT EDIT MANUALLY - Run `cargo run --package dda-spec` to regenerate.
=#

module Variants

export SPEC_VERSION, SELECT_MASK_SIZE, BINARY_NAME, REQUIRES_SHELL_WRAPPER
export SHELL_COMMAND, SUPPORTED_PLATFORMS
export ChannelFormat, Individual, Pairs, DirectedPairs
export OutputColumns, VariantMetadata
export VARIANT_REGISTRY, VARIANT_ORDER
export ST
export CT
export CD
export RESERVED
export DE
export SY

export get_variant_by_abbrev, get_variant_by_suffix, get_variant_by_position
export active_variants, generate_select_mask, parse_select_mask, format_select_mask
export ScaleParameters, generate_delays
export requires_ct_params, SelectMaskPositions
export FileType, EDF, ASCII, get_flag, file_type_from_extension

# =============================================================================
# CONSTANTS
# =============================================================================

const SPEC_VERSION = "1.0.0"
const SELECT_MASK_SIZE = 6
const BINARY_NAME = "run_DDA_AsciiEdf"
const REQUIRES_SHELL_WRAPPER = true
const SHELL_COMMAND = "sh"
const SUPPORTED_PLATFORMS = [
    "linux",
    "macos",
    "windows",
]

# =============================================================================
# ENUMS
# =============================================================================

"""Channel format for variant input."""
@enum ChannelFormat begin
    Individual
    Pairs
    DirectedPairs
end

function channel_format_from_string(s::AbstractString)::ChannelFormat
    s == "individual" && return Individual
    s == "pairs" && return Pairs
    s == "directed_pairs" && return DirectedPairs
    error("Unknown channel format: $s")
end

# =============================================================================
# TYPES
# =============================================================================

"""Output column specification for a variant."""
struct OutputColumns
    coefficients::UInt8
    has_error::Bool
    description::String
end

"""Complete variant metadata."""
struct VariantMetadata
    abbreviation::String
    name::String
    position::UInt8
    output_suffix::String
    stride::UInt8
    reserved::Bool
    required_params::Vector{String}
    channel_format::ChannelFormat
    output_columns::OutputColumns
    documentation::String
end

"""Check if variant requires CT window parameters."""
function requires_ct_params(v::VariantMetadata)::Bool
    return "-WL_CT" in v.required_params
end

# =============================================================================
# VARIANT DEFINITIONS
# =============================================================================


"""
Single Timeseries (ST) - Position 0

Analyzes individual channels independently. Most basic variant. One result row per channel.
"""
const ST = VariantMetadata(
    "ST",
    "Single Timeseries",
    0,
    "_ST",
    4,
    false,
    [],
    Individual,
    OutputColumns(
        3,
        true,
        "4 columns per channel: a_1, a_2, a_3 coefficients + error"
    ),
    "Analyzes individual channels independently. Most basic variant. One result row per channel."
)


"""
Cross-Timeseries (CT) - Position 1

Analyzes relationships between channel pairs. Symmetric: pair (1,2) equals (2,1). When enabled with ST, wrapper must run CT pairs separately.
"""
const CT = VariantMetadata(
    "CT",
    "Cross-Timeseries",
    1,
    "_CT",
    4,
    false,
    ["-WL_CT", "-WS_CT"],
    Pairs,
    OutputColumns(
        3,
        true,
        "4 columns per pair: a_1, a_2, a_3 coefficients + error"
    ),
    "Analyzes relationships between channel pairs. Symmetric: pair (1,2) equals (2,1). When enabled with ST, wrapper must run CT pairs separately."
)


"""
Cross-Dynamical (CD) - Position 2

Analyzes directed causal relationships. Asymmetric: (1->2) differs from (2->1). CD is independent (no longer requires ST+CT).
"""
const CD = VariantMetadata(
    "CD",
    "Cross-Dynamical",
    2,
    "_CD_DDA_ST",
    2,
    false,
    ["-WL_CT", "-WS_CT"],
    DirectedPairs,
    OutputColumns(
        1,
        true,
        "2 columns per directed pair: a_1 coefficient + error"
    ),
    "Analyzes directed causal relationships. Asymmetric: (1->2) differs from (2->1). CD is independent (no longer requires ST+CT)."
)


"""
Reserved (RESERVED) - Position 3

Internal development function. Should always be set to 0 in production.
"""
const RESERVED = VariantMetadata(
    "RESERVED",
    "Reserved",
    3,
    "_RESERVED",
    1,
    true,
    [],
    Individual,
    OutputColumns(
        0,
        false,
        "Reserved for internal development"
    ),
    "Internal development function. Should always be set to 0 in production."
)


"""
Delay Embedding (DE) - Position 4

Tests for ergodic behavior in dynamical systems. Produces single aggregate measure per time window (not per-channel).
"""
const DE = VariantMetadata(
    "DE",
    "Delay Embedding",
    4,
    "_DE",
    1,
    false,
    ["-WL_CT", "-WS_CT"],
    Individual,
    OutputColumns(
        0,
        false,
        "1 column: single ergodicity measure per time window"
    ),
    "Tests for ergodic behavior in dynamical systems. Produces single aggregate measure per time window (not per-channel)."
)


"""
Synchronization (SY) - Position 5

Detects synchronized behavior between signals. Produces one value per channel/measure per time window.
"""
const SY = VariantMetadata(
    "SY",
    "Synchronization",
    5,
    "_SY",
    1,
    false,
    [],
    Individual,
    OutputColumns(
        0,
        false,
        "1 column per channel/measure: synchronization coefficient"
    ),
    "Detects synchronized behavior between signals. Produces one value per channel/measure per time window."
)



"""All variants in SELECT mask order."""
const VARIANT_REGISTRY = [
    ST,
    CT,
    CD,
    RESERVED,
    DE,
    SY,
]

"""Variant abbreviations in SELECT mask order."""
const VARIANT_ORDER = [
    "ST",
    "CT",
    "CD",
    "RESERVED",
    "DE",
    "SY",
]

# =============================================================================
# SELECT MASK POSITIONS
# =============================================================================

module SelectMaskPositions
    const ST = 0
    const CT = 1
    const CD = 2
    const RESERVED = 3
    const DE = 4
    const SY = 5
end

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

"""Look up variant by abbreviation."""
function get_variant_by_abbrev(abbrev::AbstractString)::Union{VariantMetadata, Nothing}
    for v in VARIANT_REGISTRY
        v.abbreviation == abbrev && return v
    end
    return nothing
end

"""Look up variant by output suffix."""
function get_variant_by_suffix(suffix::AbstractString)::Union{VariantMetadata, Nothing}
    for v in VARIANT_REGISTRY
        v.output_suffix == suffix && return v
    end
    return nothing
end

"""Look up variant by position in SELECT mask."""
function get_variant_by_position(position::Integer)::Union{VariantMetadata, Nothing}
    for v in VARIANT_REGISTRY
        v.position == position && return v
    end
    return nothing
end

"""Get all non-reserved variants."""
function active_variants()::Vector{VariantMetadata}
    return [v for v in VARIANT_REGISTRY if !v.reserved]
end

"""
Generate a SELECT mask from variant abbreviations.

# Examples
```julia
generate_select_mask(["ST", "SY"])  # [1, 0, 0, 0, 0, 1]
```
"""
function generate_select_mask(variants::Vector{<:AbstractString})::Vector{UInt8}
    mask = zeros(UInt8, SELECT_MASK_SIZE)
    for abbrev in variants
        v = get_variant_by_abbrev(abbrev)
        if v !== nothing
            mask[v.position + 1] = 1  # Julia is 1-indexed
        end
    end
    return mask
end

"""
Parse a SELECT mask back to variant abbreviations.

# Examples
```julia
parse_select_mask([1, 0, 0, 0, 0, 1])  # ["ST", "SY"]
```
"""
function parse_select_mask(mask::Vector{<:Integer})::Vector{String}
    result = String[]
    for (pos, bit) in enumerate(mask)
        if bit == 1
            v = get_variant_by_position(pos - 1)  # Convert from 1-indexed
            if v !== nothing && !v.reserved
                push!(result, v.abbreviation)
            end
        end
    end
    return result
end

"""Format SELECT mask as space-separated string for CLI."""
function format_select_mask(mask::Vector{<:Integer})::String
    return join(string.(mask), " ")
end

# =============================================================================
# SCALE PARAMETERS
# =============================================================================

"""Parameters for generating delay values from scale range."""
@kwdef mutable struct ScaleParameters
    scale_min::Float64 = 1
    scale_max::Float64 = 20
    scale_num::Int = 20
end

"""
Generate delay values from scale parameters.

# Examples
```julia
params = ScaleParameters(scale_min=1.0, scale_max=10.0, scale_num=10)
generate_delays(params)  # [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
```
"""
function generate_delays(params::ScaleParameters)::Vector{Int}
    if params.scale_num == 1
        return [round(Int, params.scale_min)]
    end

    step = (params.scale_max - params.scale_min) / (params.scale_num - 1)
    return [round(Int, params.scale_min + i * step) for i in 0:(params.scale_num - 1)]
end

# =============================================================================
# FILE TYPES
# =============================================================================

"""File type enumeration."""
@enum FileType begin
    EDF
    ASCII
end

const FILE_TYPE_FLAGS = Dict{FileType, String}(
    EDF => "-EDF",
    ASCII => "-ASCII",
)

const FILE_TYPE_EXTENSIONS = Dict{String, FileType}(
    "edf" => EDF,
    "ascii" => ASCII,
    "txt" => ASCII,
    "csv" => ASCII,
)

"""Get CLI flag for file type."""
function get_flag(ft::FileType)::String
    return FILE_TYPE_FLAGS[ft]
end

"""Detect file type from extension."""
function file_type_from_extension(ext::AbstractString)::Union{FileType, Nothing}
    normalized = lowercase(lstrip(ext, '.'))
    return get(FILE_TYPE_EXTENSIONS, normalized, nothing)
end

end # module
