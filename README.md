# DelayDifferentialAnalysis.jl

[![CI](https://github.com/sdraeger/DelayDifferentialAnalysis.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/sdraeger/DelayDifferentialAnalysis.jl/actions/workflows/CI.yml)

Julia bindings for the DDA binary.

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

## Binary Setup

The package wraps `run_DDA_AsciiEdf`.

You can resolve the binary in two ways:

- Pass `binary_path="/full/path/to/run_DDA_AsciiEdf"` to the API call
- Rely on the existing environment/search-path fallback (`DDA_BINARY_PATH`, `~/.local/bin`, `~/bin`, `/usr/local/bin`, `/opt/dda/bin`)

## Quick Start

```julia
using DelayDifferentialAnalysis

result = run_st(
    "data.edf",
    [1, 2, 3];
    binary_path="/opt/dda/bin/run_DDA_AsciiEdf",
    model=[1, 2, 10],
    model_dimension=4,
    sampling_rate=(500, 1000),
    wl=2048,
    ws=1024,
)

println(n_channels(result))
println(n_windows(result))
println(size(result.coefficients))
```

## Generic Binary API

Use `run_analysis` when you want to call the binary directly with arbitrary variants without constructing a request object.

```julia
using DelayDifferentialAnalysis

result = run_analysis(
    "recording.edf",
    [1, 2, 3, 4],
    ["ST", "SY"];
    binary_path="/opt/dda/bin/run_DDA_AsciiEdf",
    model=[1, 2, 10],
    model_dimension=4,
    delays=[7, 10],
    window_length=2048,
    window_step=1024,
    time_range=(0.0, 50_000.0),
    sampling_rate=(500, 1000),
    out_fn=nothing,
)

println(size(result.q_matrix))
for vr in result.variant_results
    println("$(vr.variant_id): $(size(vr.q_matrix))")
end
```

The structured variant is also available:

```julia
raw = run_analysis_structured(
    "recording.edf",
    [1, 2, 3],
    ["ST", "DE"];
    binary_path="/opt/dda/bin/run_DDA_AsciiEdf",
)
```

## Variant-Specific Helpers

The package also provides typed wrappers:

- `run_st(file_path, channels; ...)`
- `run_ct(file_path, channels; ...)`
- `run_de(file_path, channels; ...)`

Use `run_ct(...)` for pairwise CT analysis across multiple channels. The generic `run_analysis(...)` function is the raw binary wrapper.

Each wrapper also accepts an in-memory `channels × samples` matrix instead of a file path:

```julia
data = randn(3, 10_000)
result = run_st(
    data;
    binary_path="/opt/dda/bin/run_DDA_AsciiEdf",
    wl=200,
    ws=100,
)
```

## Important Parameter Notes

- `channels` are 1-indexed everywhere in the Julia API
- File-based calls infer channel labels from EDF headers and from optional ASCII/TSV header rows. Pass `channel_labels` to override them explicitly.
- `model` maps directly to the binary `-MODEL` argument and defaults to `[1, 2, 10]`
- `model_dimension` is the DDA model dimension passed to `-dm`
- `sampling_rate` maps directly to `-SR low high` and defaults to `(500, 1000)`
- `out_fn` is `nothing` by default. In that case the wrapper uses a temporary output base for the call. If you pass `out_fn`, that exact value is sent to `-OUT_FN`

## Low-Level Helpers

```julia
mask = generate_select_mask(["ST", "SY"])
println(mask)                      # [1, 0, 0, 0, 0, 1]
println(format_select_mask(mask)) # "1 0 0 0 0 1"

st = get_variant_by_abbrev("ST")
println(st.name)
println(st.output_suffix)

path = find_binary("/opt/dda/bin/run_DDA_AsciiEdf")
println(path)
```

## Variants

| Abbreviation | Name              | Description                                     |
| ------------ | ----------------- | ----------------------------------------------- |
| `ST`         | Single Timeseries | Analyzes individual channels independently      |
| `CT`         | Cross-Timeseries  | Analyzes relationships between channel pairs    |
| `CD`         | Cross-Dynamical   | Analyzes directed causal relationships          |
| `DE`         | Delay Embedding   | Tests for ergodic behavior in dynamical systems |
| `SY`         | Synchronization   | Detects synchronized behavior between signals   |

## License

MIT
