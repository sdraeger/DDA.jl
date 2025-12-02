"""
    DelayDifferentialAnalysis

Julia bindings for Delay Differential Analysis (DDA).

This package provides:
- High-level `run_analysis()` function to execute DDA
- Variant metadata (ST, CT, CD, DE, SY)
- SELECT mask generation and parsing
- File type detection
- Binary discovery

# Quick Start
```julia
using DelayDifferentialAnalysis

# Run DDA analysis on an EDF file
request = DDARequest(
    "data.edf",
    [0, 1, 2],      # channels (0-based)
    ["ST", "SY"];   # variants
    window_length=2048,
    window_step=1024
)
result = run_analysis(request)

# Access results
println("Q matrix size: ", size(result.q_matrix))
for vr in result.variant_results
    println("\$(vr.variant_name): ", size(vr.q_matrix))
end
```

# Low-level API
```julia
# Generate SELECT mask for CLI
mask = generate_select_mask(["ST", "SY"])
println(format_select_mask(mask))  # "1 0 0 0 0 1"

# Look up variant metadata
st = get_variant_by_abbrev("ST")
println("ST stride: \$(st.stride)")

# Find the DDA binary
path = find_binary()
```
"""
module DelayDifferentialAnalysis

using UUIDs
using Dates

# Include the generated Variants module
include("Variants.jl")

# Re-export everything from Variants
using .Variants
export SPEC_VERSION, SELECT_MASK_SIZE, BINARY_NAME, REQUIRES_SHELL_WRAPPER
export SHELL_COMMAND, SUPPORTED_PLATFORMS
export BINARY_ENV_VAR, BINARY_HOME_ENV_VAR, DEFAULT_BINARY_PATHS
export find_binary, require_binary
export ChannelFormat, Individual, Pairs, DirectedPairs
export OutputColumns, VariantMetadata
export VARIANT_REGISTRY, VARIANT_ORDER
export ST, CT, CD, RESERVED, DE, SY
export get_variant_by_abbrev, get_variant_by_suffix, get_variant_by_position
export active_variants, generate_select_mask, parse_select_mask, format_select_mask
export DEFAULT_DELAYS
export requires_ct_params, SelectMaskPositions
export FileType, EDF, ASCII, get_flag, file_type_from_extension

# Include the Runner module
include("Runner.jl")

# Re-export Runner types and functions
using .Runner
export DDARequest, DDAResult, VariantResultData
export DDARunner, run_analysis

end # module DelayDifferentialAnalysis
