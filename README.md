# DelayDifferentialAnalysis.jl

Julia interface for the `run_DDA_AsciiEdf` binary (Cosmopolitan Libc APE format). Download it from: [https://snl.salk.edu/~claudia/DDALAB/ddalab.html](https://snl.salk.edu/~claudia/DDALAB/ddalab.html).

## Overview

This package provides a clean, type-safe Julia interface to execute the Delay Differential Analysis (DDA) binary and parse its output. It handles all the complexities of:

- Cross-platform APE binary execution (Windows/macOS/Linux)
- Command-line argument construction
- Output file parsing and matrix transformation
- Error handling and logging

## Features

- **Type-safe API**: Strongly-typed request and response structures
- **Cross-platform**: Handles APE binary execution on Unix (sh wrapper) and Windows (direct exe)
- **Automatic parsing**: Processes raw DDA output into usable matrices
- **Error handling**: Comprehensive error types with descriptive messages

## Installation

```julia
using Pkg
Pkg.add("DelayDifferentialAnalysis")
```

For development version:

```julia
using Pkg
Pkg.add(url="https://github.com/sdraeger/DelayDifferentialAnalysis.jl")
```

## Usage

```julia
using DelayDifferentialAnalysis

# Create runner with path to run_DDA_AsciiEdf binary
runner = DDARunner("/path/to/run_DDA_AsciiEdf")

# Build analysis request
request = DDARequest(
    "/path/to/data.edf",                     # file_path
    [1, 2],                                  # channels (1-based indices)
    nothing,                                 # bounds (use nothing to process entire file, or Bounds(start, stop) for a range)
    AlgorithmSelection(
        ["standard"],                         # enabled_variants
        nothing                               # select_mask
    ),
    WindowParameters(
        1024,                                 # window_length
        512,                                  # window_step
        nothing,                              # ct_window_length
        nothing                               # ct_window_step
    ),
    DelayParameters([1, 2, 3, 4, 5, 6, 7, 8, 9, 10]),  # delays (can be any list of integers)
    nothing                                   # ct_channel_pairs
)

# Run analysis
result = run_dda(runner, request)

# Access results
println("Q matrix shape: $(size(result.q_matrix))")
```

## Architecture

### Components

- **`types.jl`**: Request/response structures and parameter types
- **`runner.jl`**: DDA binary execution logic
- **`parser.jl`**: Output file parsing and matrix transformation
- **`error.jl`**: Error types and custom exceptions

### Binary Execution

The package automatically handles the APE (Actually Portable Executable) format:

- **Unix (macOS/Linux)**: Runs through `sh` wrapper to handle polyglot format
- **Windows**: Executes `.exe` directly

### Output Processing

The parser implements the same transformation as dda-py and dda-rs:

1. Skip first 2 columns
2. Take every 4th column from the remaining data
3. Transpose to get [channels × timepoints] format

## API Reference

### Types

#### `DDARunner`

Main runner type that wraps the binary path.

#### `DDARequest`

Complete configuration for DDA analysis.

#### `DDAResult`

Analysis results containing Q matrices and metadata.

#### `AlgorithmSelection`

Specifies which DDA variants to enable. Contains:
- `enabled_variants`: Vector of variant names (e.g., `["ST", "CT", "CD", "DE"]`)
- `select_mask`: Optional manual select mask string (e.g., `"1 1 0 0"`). If `nothing`, the mask is auto-generated from `enabled_variants`.

**Note**: CD variant automatically enables ST and CT as they are required dependencies.

#### `Bounds`

Sample bounds specification with `start` and `stop` fields (integers representing sample indices). Use `nothing` to process the entire file, or provide `Bounds(start, stop)` to process a specific sample range.

#### `WindowParameters`

Window configuration for DDA analysis.

#### `DelayParameters`

Delay parameter configuration. Accepts a vector of integer delay values (e.g., `[1, 2, 3, 5, 10]` or `collect(1:20)`).

### Functions

#### `DDARunner(binary_path::String)`

Create a new DDA runner with the specified binary path.

**Throws:** `BinaryNotFoundError` if the binary doesn't exist.

#### `run_dda(runner::DDARunner, request::DDARequest, edf_channel_names::Union{Vector{String},Nothing}=nothing)`

Run DDA analysis with the given request parameters. Sample bounds are specified within the `request.bounds` field (use `nothing` to process the entire file).

**Returns:** `DDAResult` containing the processed Q matrix and metadata.

**Throws:**

- `FileNotFoundError` if input file doesn't exist
- `ExecutionFailedError` if binary execution fails
- `ParseError` if output parsing fails

#### `binary_path(runner::DDARunner)`

Get the path to the DDA binary.

### Error Types

All errors inherit from `DDAError`:

- `BinaryNotFoundError`: DDA binary not found
- `FileNotFoundError`: Input file not found
- `UnsupportedFileTypeError`: Unsupported file type
- `ExecutionFailedError`: Binary execution failed
- `ParseError`: Failed to parse output
- `InvalidParameterError`: Invalid parameter provided

## Examples

### Basic Single-Channel Analysis

```julia
using DelayDifferentialAnalysis

runner = DDARunner("./bin/run_DDA_AsciiEdf")

request = DDARequest(
    "data.edf",
    [1, 2, 3],  # First 3 channels (1-based indices)
    nothing,    # Process entire file
    AlgorithmSelection(["standard"], nothing),
    WindowParameters(512, 256, nothing, nothing),
    DelayParameters(collect(1:20)),  # Or specify custom delays like [1, 5, 10, 15, 20]
    nothing
)

result = run_dda(runner, request)
println("Q matrix: $(size(result.q_matrix))")
```

### Multi-Channel Analysis with CT Variant and Sample Bounds

```julia
using DelayDifferentialAnalysis

runner = DDARunner("./bin/run_DDA_AsciiEdf")

request = DDARequest(
    "data.edf",
    [1, 2, 3],
    Bounds(0, 60000),  # Process samples 0-60000
    AlgorithmSelection(["ST", "CT"], nothing),  # Enable ST and CT (select_mask auto-generated)
    WindowParameters(1024, 512, 2048, 1024),  # CT-specific windows
    DelayParameters(collect(1:15)),  # Delays from 1 to 15
    [[1, 2], [2, 3]]  # CT channel pairs
)

result = run_dda(runner, request)

# Access variant-specific results
if !isnothing(result.variant_results)
    for variant in result.variant_results
        println("$(variant.variant_name): $(size(variant.q_matrix))")
    end
end
```

### Cross-Dynamical (CD) Analysis

```julia
using DelayDifferentialAnalysis

runner = DDARunner("./bin/run_DDA_AsciiEdf")

request = DDARequest(
    "data.edf",
    [1, 2],
    nothing,  # Process entire file
    AlgorithmSelection(["CD"], nothing),  # Enable CD (auto-enables ST and CT too)
    WindowParameters(1024, 512, nothing, nothing),
    DelayParameters([1, 2, 5, 10]),  # Custom delay values
    nothing
)

result = run_dda(runner, request)

# CD variant uses a different output format with stride=2 for column extraction
if !isnothing(result.variant_results)
    for variant in result.variant_results
        println("$(variant.variant_name): $(size(variant.q_matrix))")
    end
end
```

## Testing

```julia
using Pkg
Pkg.test("DelayDifferentialAnalysis")
```

## License

MIT

## Contributing

Contributions are welcome! Please ensure all tests pass before submitting a pull request.

## Related Packages

- **dda-rs**: Rust interface for run_DDA_AsciiEdf
- **dda-py**: Python interface for run_DDA_AsciiEdf
