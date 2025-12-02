"""
    DelayDifferentialAnalysis

Julia bindings for Delay Differential Analysis (DDA).

This package provides:
- Variant metadata (ST, CT, CD, DE, SY)
- SELECT mask generation and parsing
- File type detection
- Scale parameter utilities

# Example
```julia
using DelayDifferentialAnalysis

# Generate SELECT mask for ST and SY variants
mask = generate_select_mask(["ST", "SY"])

# Look up variant metadata
st = get_variant_by_abbrev("ST")
println("ST stride: \$(st.stride)")
```
"""
module DelayDifferentialAnalysis

# Include the generated Variants module
include("Variants.jl")

# Re-export everything from Variants
using .Variants
export SPEC_VERSION, SELECT_MASK_SIZE, BINARY_NAME, REQUIRES_SHELL_WRAPPER
export SHELL_COMMAND, SUPPORTED_PLATFORMS
export ChannelFormat, Individual, Pairs, DirectedPairs
export OutputColumns, VariantMetadata
export VARIANT_REGISTRY, VARIANT_ORDER
export ST, CT, CD, RESERVED, DE, SY
export get_variant_by_abbrev, get_variant_by_suffix, get_variant_by_position
export active_variants, generate_select_mask, parse_select_mask, format_select_mask
export ScaleParameters, generate_delays
export requires_ct_params, SelectMaskPositions
export FileType, EDF, ASCII, get_flag, file_type_from_extension

end # module DelayDifferentialAnalysis
