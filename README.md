# DDA.jl

Julia interface for the `run_DDA_AsciiEdf` binary (Cosmopolitan Libc APE format).

## Overview

This package provides a clean, type-safe Julia interface to execute the DDA (Delay Differential Analysis) binary and parse its output. It handles all the complexities of:

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
Pkg.add("DDA")
```

For development version:

```julia
using Pkg
Pkg.add(url="https://github.com/USERNAME/DDA.jl")
```

## Usage

```julia
using DDA

# Create runner with path to run_DDA_AsciiEdf binary
runner = DDARunner("/path/to/run_DDA_AsciiEdf")

# Build analysis request
request = DDARequest(
    "/path/to/data.edf",                    # file_path
    [0, 1],                                  # channels (0-based indices)
    TimeRange(0.0, 100.0),                   # time_range
    PreprocessingOptions(
        "linear",                             # detrending
        nothing,                              # highpass
        nothing                               # lowpass
    ),
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
    ScaleParameters(
        1.0,                                  # scale_min
        10.0,                                 # scale_max
        10                                    # scale_num
    ),
    nothing                                   # ct_channel_pairs
)

# Run analysis with sample bounds (start_sample, end_sample)
start_bound = 0        # Start from beginning
end_bound = 10000      # End at sample 10000
result = run(runner, request, start_bound, end_bound)

# Access results
println("Analysis ID: $(result.id)")
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
3. Transpose to get [channels/scales × timepoints] format

## API Reference

### Types

#### `DDARunner`
Main runner type that wraps the binary path.

#### `DDARequest`
Complete configuration for DDA analysis.

#### `DDAResult`
Analysis results containing Q matrices and metadata.

#### `TimeRange`
Time range specification with `start` and `stop` fields.

#### `WindowParameters`
Window configuration for DDA analysis.

#### `ScaleParameters`
Scale/delay parameter configuration.

### Functions

#### `DDARunner(binary_path::String)`
Create a new DDA runner with the specified binary path.

**Throws:** `BinaryNotFoundError` if the binary doesn't exist.

#### `run(runner::DDARunner, request::DDARequest, start_bound::Int, end_bound::Int, edf_channel_names::Union{Vector{String},Nothing}=nothing)`
Run DDA analysis with the given request parameters.

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
using DDA

runner = DDARunner("./bin/run_DDA_AsciiEdf")

request = DDARequest(
    "data.edf",
    [0],  # First channel
    TimeRange(0.0, 60.0),
    PreprocessingOptions("linear", nothing, nothing),
    AlgorithmSelection(["standard"], nothing),
    WindowParameters(512, 256, nothing, nothing),
    ScaleParameters(1.0, 20.0, 20),
    nothing
)

result = run(runner, request, 0, 30000)
println("Q matrix: $(size(result.q_matrix))")
```

### Multi-Channel Analysis with CT Variant

```julia
using DDA

runner = DDARunner("./bin/run_DDA_AsciiEdf")

request = DDARequest(
    "data.edf",
    [0, 1, 2],
    TimeRange(0.0, 120.0),
    PreprocessingOptions("linear", 0.5, 45.0),
    AlgorithmSelection(["ST", "CT"], "1 1 0 0"),  # Enable ST and CT
    WindowParameters(1024, 512, 2048, 1024),  # CT-specific windows
    ScaleParameters(1.0, 15.0, 15),
    [[0, 1], [1, 2]]  # CT channel pairs
)

result = run(runner, request, 0, 60000)

# Access variant-specific results
if !isnothing(result.variant_results)
    for variant in result.variant_results
        println("$(variant.variant_name): $(size(variant.q_matrix))")
    end
end
```

## Testing

```julia
using Pkg
Pkg.test("DDA")
```

## License

MIT

## Contributing

Contributions are welcome! Please ensure all tests pass before submitting a pull request.

## Related Packages

- **dda-rs**: Rust interface for run_DDA_AsciiEdf
- **dda-py**: Python interface for run_DDA_AsciiEdf
