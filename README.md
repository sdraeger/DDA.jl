# DelayDifferentialAnalysis.jl

[![CI](https://github.com/sdraeger/DelayDifferentialAnalysis.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/sdraeger/DelayDifferentialAnalysis.jl/actions/workflows/CI.yml)

Julia package for Delay Differential Analysis (DDA).

## Requirements

The [DDA binary](https://snl.salk.edu/~sfdraeger/dda/) (`run_DDA_AsciiEdf`) is required. Download the most recent version from the file server and place it in one of these locations:

- Set `$DDA_BINARY_PATH` environment variable to the full path
- Set `$DDA_HOME` and place binary in `$DDA_HOME/bin/`
- Place in `~/.local/bin/`, `~/bin/`, `/usr/local/bin/`, or `/opt/dda/bin/`

## Installation

```julia
using Pkg
Pkg.add("DelayDifferentialAnalysis")
```

Or, for the development version:

```julia
using Pkg
Pkg.add(url="https://github.com/sdraeger/DelayDifferentialAnalysis.jl")
```

## Quick Start

```julia
using DelayDifferentialAnalysis

# Create analysis request
request = DDARequest(
    "data.edf",           # Input file (EDF or ASCII)
    [0, 1, 2],            # Channel indices (0-based)
    ["ST", "SY"];         # Variants to run
    window_length=2048,
    window_step=1024
)

# Run analysis (auto-discovers binary)
result = run_analysis(request)

# Access results
println("Q matrix size: ", size(result.q_matrix))

for vr in result.variant_results
    println("$(vr.variant_name): $(size(vr.q_matrix))")
end
```

## Full Example

```julia
using DelayDifferentialAnalysis

# Create request with all parameters
request = DDARequest(
    "recording.edf",
    [0, 1, 2, 3],              # Channels (0-based)
    ["ST", "CT", "SY"];        # Variants
    window_length=2048,        # Analysis window length
    window_step=1024,          # Window step size
    ct_window_length=4,        # CT-specific window length
    ct_window_step=2,          # CT-specific window step
    delays=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10],  # Delay (tau) values
    dm=4,                      # Embedding dimension
    order=4,                   # Polynomial order
    nr_tau=2,                  # Number of tau values
    time_range=(0.0, 50000.0), # Sample range (optional)
    ct_pairs=[(0, 1), (0, 2)], # CT channel pairs (optional)
    sampling_rate=256.0        # Sampling rate in Hz (optional)
)

# Option 1: Auto-discover binary
result = run_analysis(request)

# Option 2: Explicit binary path
runner = DDARunner("/path/to/run_DDA_AsciiEdf")
result = run_analysis(runner, request)

# Access variant-specific results
for vr in result.variant_results
    println("Variant: $(vr.variant_id)")
    println("  Name: $(vr.variant_name)")
    println("  Q matrix: $(size(vr.q_matrix))")
    println("  Channels: $(vr.channel_labels)")
end
```

## Low-Level API

```julia
using DelayDifferentialAnalysis

# Generate SELECT mask for DDA CLI
mask = generate_select_mask(["ST", "SY"])
println(mask)                    # [1, 0, 0, 0, 0, 1]
println(format_select_mask(mask)) # "1 0 0 0 0 1"

# Parse SELECT mask back to variant names
variants = parse_select_mask([1, 1, 0, 0, 0, 1])
println(variants)  # ["ST", "CT", "SY"]

# Look up variant metadata
st = get_variant_by_abbrev("ST")
println("Name: $(st.name)")           # "Single Timeseries"
println("Stride: $(st.stride)")       # 4
println("Suffix: $(st.output_suffix)") # "_ST"

# Check if variant requires CT parameters
println(requires_ct_params(CT))  # true
println(requires_ct_params(ST))  # false

# File type detection
ft = file_type_from_extension(".edf")
println(get_flag(ft))  # "-EDF"

# Find binary location
path = find_binary()
if path !== nothing
    println("Binary found at: $path")
end
```

## Variants

| Abbreviation | Name | Description |
|--------------|------|-------------|
| `ST` | Single Timeseries | Analyzes individual channels independently |
| `CT` | Cross-Timeseries | Analyzes relationships between channel pairs |
| `CD` | Cross-Dynamical | Analyzes directed causal relationships |
| `DE` | Delay Embedding | Tests for ergodic behavior in dynamical systems |
| `SY` | Synchronization | Detects synchronized behavior between signals |

## API Reference

### Types

- `DDARequest` - Analysis configuration
- `DDAResult` - Analysis results with Q matrices
- `VariantResultData` - Per-variant result data
- `DDARunner` - Binary wrapper

### Functions

- `run_analysis(request)` - Run DDA with auto-discovered binary
- `run_analysis(runner, request)` - Run DDA with explicit runner
- `DDARunner()` - Create runner with auto-discovery
- `DDARunner(path)` - Create runner with explicit path
- `find_binary()` - Find DDA binary location
- `require_binary()` - Find binary or throw error
- `generate_select_mask(variants)` - Create SELECT mask
- `parse_select_mask(mask)` - Parse mask to variant names
- `format_select_mask(mask)` - Format mask for CLI
- `get_variant_by_abbrev(abbrev)` - Look up variant by abbreviation
- `active_variants()` - Get all non-reserved variants

## License

MIT
