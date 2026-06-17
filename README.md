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

On Windows, rename the binary to include the `.exe` suffix and pass that path,
for example `binary_path="C:\\path\\to\\run_DDA_AsciiEdf.exe"`.

## Quick Start

```julia
using DelayDifferentialAnalysis

result = run_st(
    file_path="data.edf",
    channels=[1, 2, 3],
    binary_path="/opt/dda/bin/run_DDA_AsciiEdf",
    derivative_points=3,
    WL=2048,
    WS=1024,
)

println(n_channels(result))
println(n_windows(result))
println(size(result.coefficients))
println(result.T[1:3])
println(result.t[1:3])
```

## Generic Binary API

Use `run_DDA` when you want to call the binary directly with arbitrary flavors without constructing a request object.

```julia
using DelayDifferentialAnalysis

result = run_DDA(
    file_path="recording.edf",
    channels=[1, 2, 3, 4],
    flavors=["ST", "SY"],
    binary_path="/opt/dda/bin/run_DDA_AsciiEdf",
    model=[1, 2, 10],
    derivative_points=3,
    order=4,
    delays=[7, 10],
    WL=2048,
    WS=1024,
    time_range=(0, 50000),
    out_fn=nothing,
)

println(size(result.T))  # first two raw integer binary columns
println(size(result.t))  # derived time axis from result.T[:, 1]
println(size(result.A))  # all remaining raw binary columns
println(size(result.ST))
println(size(result.SY))
```

`channels` is optional for `run_DDA`. If it is omitted or set to `nothing`, the
wrapper does not pass `-CH_list` and leaves channel selection to the binary.
Pass `load_results=false` to execute the binary without parsing output files
into Julia; the call returns `nothing`. If output files should be kept in this
mode, pass an explicit `out_fn`.

## Structure Selection

`make_MOD` generates the DDA model library `MOD` used by the original
structure-selection scripts. `N_MOD` is the number of monomial terms per model,
`DDAorder` is the polynomial order, and `nr_delays` defaults to `2`.

```julia
MOD = make_MOD(3, 3)
print_structure_selection(MOD, 3)
```

The printer writes the `P_DDA` monomial encoding table, a checkmark table where
each `MOD` column corresponds to a `P_DDA` monomial, and each `MOD` row as a
Unicode equation for direct terminal inspection, for example
`ẋ = a₁·x₁ + a₂·x₂²`. The low-level terminal printer is available as
`write_model_terminal(...)`; `write_model_LaTeX(...)` remains available when
LaTeX output is needed.

`structure_selection` evaluates candidate `-MODEL` and `-TAU` combinations
with the existing `run_DDA` path and returns the lowest-error candidate. Pass
explicit `candidate_models`, an existing `MOD`, or `N_MOD` plus `DDAorder` to
generate candidates with `make_MOD`. Pass `delays` as a flat delay pool to
write Claudia-style `TAU_ALL__...` files for each model symmetry class, or pass
nested delay vectors for explicit delay candidates.

```julia
selection = structure_selection(
    file_path="recording.ascii",
    channels=[1, 2],
    binary_path="/opt/dda/bin/run_DDA_AsciiEdf",
    N_MOD=3,
    DDAorder=3,
    delays=(4 + 1):40,
    derivative_points=4,
    WL=3000,
    WS=200,
    input_format=:ascii,
    tau_file_prefix="/scratch/my-run/TAU_ALL__",
)

println(selection.best_model)
println(selection.best_delays)
println(selection.best_score)
```

By default, `structure_selection` uses `model_scope=:joint`, which selects one
model/delay combination for all requested channels together. Pass
`model_scope=:per_channel` or `model_scope="per_channel"` to run the same
candidate search independently for each channel and return one selected model
per channel.

Candidate models use the same encoding as `run_DDA`: either vectors of binary
`-MODEL` indices or matrices whose rows are monomial encodings. The old
`candidate_delays` keyword remains as an alias for nested explicit delay
candidates, but new code should use `delays`. A flat vector or range such as
`delays=[10, 20, 30]` or `delays=(derivative_points + 1):TM` is treated as a
delay pool and written to model-specific tau files in one generated artifact
folder per structure-selection call. Pass `out_dir` to control where that
folder is created, or pass `tau_file_prefix` to choose the exact path prefix
for tau files. For example, `tau_file_prefix="/scratch/my-run/TAU_ALL__"`
writes files such as `/scratch/my-run/TAU_ALL__1_0` and
`/scratch/my-run/TAU_ALL__2_1`, avoiding systems where `/tmp` is not writable.
Artifacts are retained by default for external inspection; pass
`cleanup_on_error=true` to delete only automatically generated artifact folders
if the selection call errors.

The structured variant is also available:

```julia
raw = run_analysis_structured(
    file_path="recording.edf",
    channels=[1, 2, 3],
    flavors=["ST", "DE"],
    binary_path="/opt/dda/bin/run_DDA_AsciiEdf",
)
```

## Variant-Specific Helpers

The package also provides typed wrappers:

- `run_st(file_path=..., channels=...; ...)`
- `run_ct(file_path=..., channels=...; ...)`
- `run_de(file_path=..., channels=...; ...)`

Use `run_ct(...)` for pairwise CT analysis across multiple channels. The generic `run_DDA(...)` function is the raw binary wrapper.

Each wrapper also accepts an in-memory `channels × samples` matrix instead of a file path:

```julia
data = randn(3, 10_000)
result = run_st(
    data=data;
    binary_path="/opt/dda/bin/run_DDA_AsciiEdf",
    WL=200,
    WS=100,
)
```

## Important Parameter Notes

- `channels` are 1-indexed everywhere in the Julia API
- File-based calls infer channel labels from EDF headers and from optional ASCII/TSV header rows. Pass `channel_labels` to override them explicitly.
- File-based calls infer the binary input flag from the extension by default: `.edf` emits `-EDF`, all other extensions emit `-ASCII`. Pass `input_format=:ascii` or `input_format=:edf` to override this manually; string values such as `input_format="ascii"` are also accepted.
- `model` maps to the binary `-MODEL` argument. Pass a vector of `-MODEL` indices directly, or pass a matrix whose rows are monomial encodings; matrix rows are converted to indices using `nr_tau` and `order`. If you pass a custom model, also pass explicit `derivative_points` and `order`
- `derivative_points` is the preferred Julia name for the binary `-dm` parameter and defaults to `3`
- `WL` and `WS` map to the binary `-WL` and `-WS` arguments. Both default to `nothing`; unset values are not passed to the binary.
- `-WLms` and `-WSms` are special binary flags and are intentionally not emitted by this wrapper.
- `WL_CT` and `WS_CT` are channel-group parameters, not temporal window aliases. `run_DDA` emits each flag at most once. The legacy aliases `ct_window_length` and `ct_window_step` are accepted when they agree with `WL_CT` and `WS_CT`, and conflicting values raise an error.
- `run_DDA` accepts raw passthrough keywords for advanced binary options: `tau_file::String` maps to `-TAU_file`, `tau2::Vector{Int}` maps to `-TAU2`, `model2::Vector{Int}` maps to `-MODEL2`, `no_norm::Bool` maps to `-NoNorm`, and `WN_list::Vector{Int}` maps to `-WN_list`. These default to `nothing` or `false` and are not passed unless specified. When `tau_file` is specified, `delays` are ignored for shell command generation and no direct `-TAU` list is passed.
- `run_DDA` executes the requested binary command once and parses the native output files produced by the binary, including mixed `ST`/`CT` runs.
- `select` can be passed to `run_DDA(...)` or `run_analysis_structured(...)` as a raw `-SELECT` mask. When present, it overrides the string `flavors` list
- Low-level `run_DDA` results expose each returned flavor matrix directly as a property, for example `result.ST` and `result.CT`. Top-level `result.T`, `result.t`, and `result.A` mirror the first returned flavor for backward compatibility.
- `TM` is used only for `result.t` and defaults to `max(delays)`
- `sampling_rate` is optional and defaults to `nothing`. When omitted, no `-SR` flag is passed. A scalar maps to `-SR N`; a tuple maps to `-SR N1 N2`.
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

## Flavors

| Abbreviation | Name                 | Description                                     
| ------------ | -----------------    | ----------------------------------------------- 
| `ST`         | Single Timeseries    | Analyzes individual channels independently      
| `CT`         | Cross-Timeseries     | DDA on multiple channels simultaniously   
| `CD`         | Cross-Dynamical      | Analyzes directed causal relationships          
| `DE`         | Dynamical Ergodicity | Tests for dynamical similarity in dynamical systems 
| `SY`         | Synchronization      | Detects synchronized behavior between signals   

## License

MIT
