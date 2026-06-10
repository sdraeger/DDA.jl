"""
    DelayDifferentialAnalysis

Julia bindings for Delay Differential Analysis (DDA).

# High-Level API
```julia
using DelayDifferentialAnalysis

result = run_st(
    file_path="data.edf",
    channels=[1, 2, 3],
    binary_path="/opt/dda/bin/run_DDA_AsciiEdf",
    derivative_points=3,
    WL=200,
    WS=100,
)
println(n_channels(result))           # 3
println(n_windows(result))            # depends on data length
println(result.coefficients |> size)  # (3, n_windows, 3)
println(result.T[1])                  # raw first output column
println(result.t[1])                  # derived time axis
```

# Cross-Timeseries
```julia
result = run_ct(
    file_path="data.edf",
    channels=[1, 2, 3],
    binary_path="/opt/dda/bin/run_DDA_AsciiEdf",
    WL=200,
    WS=100,
)
println(n_pairs(result))  # 3 pairs for 3 channels
```

# Dynamical Ergodicity
```julia
result = run_de(
    file_path="data.edf",
    channels=[1, 2, 3],
    binary_path="/opt/dda/bin/run_DDA_AsciiEdf",
    WL=200,
    WS=100,
)
println(result.ergodicity |> length)
```

# Model Encoding
```julia
monomials = generate_monomials(2, 4)
println(decode_model_encoding([1, 2, 10]; num_delays=2, polynomial_order=4))
```

# Low-Level API
```julia
result = run_DDA(
    file_path="data.edf",
    channels=[1, 2, 3],
    flavors=["ST"],
    binary_path="/opt/dda/bin/run_DDA_AsciiEdf",
    derivative_points=3,
    WL=200,
    WS=100,
)
```
"""
module DelayDifferentialAnalysis

# --- Flavors (binary metadata, SELECT masks, file types) ---
include("Flavors.jl")
using .Flavors
const Variants = Flavors
export SPEC_VERSION, SELECT_MASK_SIZE, BINARY_NAME, REQUIRES_SHELL_WRAPPER
export SHELL_COMMAND, SUPPORTED_PLATFORMS
export BINARY_ENV_VAR, DEFAULT_BINARY_PATHS
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

# --- Defaults (parameter constants) ---
include("Defaults.jl")
using .DDADefaults
using .DDAFlags

# --- Model encoding ---
include("ModelEncoding.jl")
using .ModelEncoding
export generate_monomials, model_matrix_to_encoding, monomial_to_text, monomial_to_latex
export decode_model_encoding, visualize_model_space

# --- Result types ---
include("Results.jl")
using .Results
export STResult, CTResult, DEResult
export n_channels, n_windows, n_coeffs, n_pairs, to_dataframe

# --- Runner (binary execution + parsing) ---
include("Runner.jl")
using .Runner
export DDARequest, DDAResult, VariantResultData
export DDARunner, run_DDA
export StructuredTimepoint, StructuredChannelData
export run_analysis_structured, parse_output_file_structured

# --- Structure selection ---
include("StructureSelection.jl")
using .StructureSelection
export StructureSelectionTrial, StructureSelectionResult, structure_selection

# --- High-level API ---
include("API.jl")
using .API
export run_st, run_ct, run_de

# --- Batch processing ---
include("Batch.jl")
using .Batch
export GroupResult, run_batch, collect_results
export n_subjects, mean_over_windows

# --- Statistics ---
include("Stats.jl")
using .Stats
export PermutationResult, EffectSizeResult, WindowComparisonResult
export permutation_test, compute_effect_size, compare_windows

# --- Plotting (lazy Plots.jl) ---
include("Plotting.jl")
using .Plotting
export plot_coefficients, plot_heatmap, plot_errors, plot_ergodicity, plot_model

end # module DelayDifferentialAnalysis
