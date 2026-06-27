# DelayDifferentialAnalysis.jl

[![CI](https://github.com/sdraeger/DelayDifferentialAnalysis.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/sdraeger/DelayDifferentialAnalysis.jl/actions/workflows/CI.yml)

Julia bindings for the native DDA binary `run_DDA_AsciiEdf`.

The package is intentionally small: it builds the binary command, runs it, and
parses the native output files into Julia objects when requested. It does not
reimplement the DDA algorithm.

## Installation

```julia
using Pkg
Pkg.add("DelayDifferentialAnalysis")
```

For the development version:

```julia
using Pkg
Pkg.add(url="https://github.com/sdraeger/DelayDifferentialAnalysis.jl")
```

## Binary Setup

Pass the binary explicitly when possible:

```julia
binary_path = "/opt/dda/bin/run_DDA_AsciiEdf"
```

If `binary_path` is omitted, the package checks `DDA_BINARY_PATH` and common
locations such as `~/.local/bin`, `~/bin`, `/usr/local/bin`, and `/opt/dda/bin`.

On Windows, rename the binary with an `.exe` suffix and pass that path, for
example `binary_path="C:\\path\\to\\run_DDA_AsciiEdf.exe"`.

## Running DDA

Use `run_DDA` for direct access to the binary. Arguments are keyword-only.

```julia
using DelayDifferentialAnalysis

result = run_DDA(
    file_path="recording.edf",
    channels=[1, 2, 3, 4],
    flavors=["ST", "CT"],
    binary_path="/opt/dda/bin/run_DDA_AsciiEdf",
    model=[1, 2, 10],
    derivative_points=4,
    order=3,
    delays=[7, 10],
    WL=2048,
    WS=1024,
    time_range=(0, 50_000),
)

println(size(result.T))   # first two integer columns from the native output
println(size(result.t))   # derived time axis
println(size(result.A))   # remaining coefficient/error columns
println(size(result.ST))  # flavor-specific output matrix
println(size(result.CT))
```

`channels` is optional. If omitted or set to `nothing`, no `-CH_list` is passed
and channel selection is left to the binary.

Use `load_results=false` when the binary should run without loading output files
into Julia:

```julia
run_DDA(
    file_path="recording.ascii",
    flavors=["ST"],
    binary_path="/opt/dda/bin/run_DDA_AsciiEdf",
    input_format=:ascii,
    model=[1, 2, 10],
    derivative_points=4,
    order=3,
    delays=[7, 10],
    out_fn="dda-output/run1",
    load_results=false,
)
```

The typed helpers `run_st`, `run_ct`, and `run_de` call the same execution path
with common flavor defaults.

## Structure Selection

Structure selection has two recommended steps:

1. `structure_selection_compute` computes and caches DDA outputs.
2. `structure_selection_select` reads cached outputs and chooses the lowest-error
   model/delay combination.

This split is useful on shared machines because expensive binary runs can be
reused, inspected, and selected over different channel subsets without rerunning
DDA.

### Build and Inspect MOD

`make_MOD(N_MOD, DDAorder; nr_delays=2)` creates Claudia-style binary model
matrices. Each row is one candidate model. Each column corresponds to one
monomial in the printed `P_DDA` table.

```julia
MOD = make_MOD(3, 3)
print_structure_selection(MOD, 3)
```

The terminal printer shows `P_DDA`, a compact checkmark table for `MOD`, and a
Unicode model equation for each row. `write_model_LaTeX` remains available when
LaTeX output is needed explicitly.

### Compute Cached DDA Outputs

```julia
run = structure_selection_compute(
    file_path="recording.ascii",
    channels=collect(1:11),
    binary_path="/opt/dda/bin/run_DDA_AsciiEdf",
    input_format=:ascii,
    N_MOD=3,
    DDAorder=3,
    delays=(4 + 1):40,
    derivative_points=4,
    WL=3000,
    WS=200,
    prefix="RUN1",
)
```

If `MOD` is omitted, `N_MOD` and `DDAorder` are used to call `make_MOD`. If
`MOD_numbers` is omitted, every row of `MOD` is computed. If `channels` is
omitted, no `-CH_list` is passed and the binary runs over its default channel
set.

`prefix` is a folder. It is created if needed. The compute step writes generated
tau files and DDA outputs there, for example:

```text
RUN1/TAU_ALL__1_0
RUN1/TAU_ALL__2_1
RUN1/structure_selection_01_02_10_ST
RUN1/structure_selection_01_02_10.info
```

Output names encode the active `MOD` column indices as two-digit values joined
by underscores. Existing complete outputs are reused rather than overwritten.
Candidates are evaluated in randomized order by default; pass `randomize=false`
for deterministic order.

### Select From Cached Outputs

```julia
selection = structure_selection_select(run; channels=[1])

println(selection.best_model)
println(selection.best_delays)
println(selection.best_score)
```

`structure_selection_select` does not run the DDA binary. It only reads files in
`run.prefix`. If `channels` is omitted, all channels found in the cached output
are scored together. If a compute run used explicit channel IDs, selection maps
those same IDs back to their cached positions. `channel=[...]` is accepted as a
singular alias for `channels=[...]`.

Pass `model_scope=:per_channel` to choose one model/delay combination per
channel:

```julia
per_channel = structure_selection_select(run; model_scope=:per_channel)
```

The old one-shot `structure_selection(...)` wrapper is still available for
small interactive runs. New long-running workflows should prefer
`structure_selection_compute` plus `structure_selection_select`.

## Parameter Notes

- Julia channel indices are 1-based.
- `input_format` accepts `:ascii`, `:edf`, `"ascii"`, or `"edf"`. If omitted,
  `.edf` files emit `-EDF`; other extensions emit `-ASCII`.
- `model` maps to `-MODEL`. A vector is passed as model indices. A matrix is
  interpreted row-wise as monomial encodings and converted to model indices.
- `derivative_points` is the Julia name for the binary `-dm` parameter.
- `WL` and `WS` map to `-WL` and `-WS`. They default to `nothing` and are not
  passed unless specified.
- `WL_CT` and `WS_CT` map to `-WL_CT` and `-WS_CT`.
- `tau_file` maps to `-TAU_file`. When provided, direct `delays` are ignored
  for command generation and no `-TAU` list is passed.
- `tau2`, `model2`, `no_norm`, and `WN_list` pass through to `-TAU2`,
  `-MODEL2`, `-NoNorm`, and `-WN_list` when specified.
- `sampling_rate=nothing` is the default and does not pass `-SR`. A scalar maps
  to `-SR N`; a tuple maps to `-SR N1 N2`.
- `select=[...]` passes a raw `-SELECT` mask and overrides `flavors`.
- `out_fn=nothing` uses a temporary output base. Pass `out_fn` to keep output
  files at a specific path.

## Result Objects

Parsed `run_DDA` results expose:

- `result.T`: first two integer columns of the native output.
- `result.t`: derived time axis, computed from `T[:, 1]`, `TM`,
  `derivative_points`, and `sampling_rate` when available.
- `result.A`: remaining native output columns for the first returned flavor.
- `result.ST`, `result.CT`, `result.CD`, `result.DE`, `result.SY`: direct
  flavor matrices when those outputs are present.

File-based calls infer labels from EDF headers or optional ASCII/TSV header
rows. Pass `channel_labels` to override labels explicitly.

## Low-Level Helpers

```julia
mask = generate_select_mask(["ST", "SY"])
println(format_select_mask(mask))  # "1 0 0 0 0 1"

st = get_variant_by_abbrev("ST")
println(st.name)

path = find_binary("/opt/dda/bin/run_DDA_AsciiEdf")
println(path)
```

## Flavors

| Abbreviation | Name                 |
| ------------ | -------------------- |
| `ST`         | Single Timeseries    |
| `CT`         | Cross-Timeseries     |
| `CD`         | Cross-Dynamical      |
| `DE`         | Dynamical Ergodicity |
| `SY`         | Synchronization      |

## License

MIT
