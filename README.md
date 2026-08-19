# DelayDifferentialAnalysis.jl

[![CI](https://github.com/sdraeger/DelayDifferentialAnalysis.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/sdraeger/DelayDifferentialAnalysis.jl/actions/workflows/CI.yml)

Julia interfaces for Delay Differential Analysis (DDA).

The package supports both the established `run_DDA_AsciiEdf` binary and a native
Julia implementation of the ST, CT, CD, DE, and SY flavors. The native engine
can run on the CPU or use CUDA for its batched regression work.

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

## Native Julia Engine

`run_dda_matrix` runs DDA directly on a `samples × channels` matrix without an
external binary. Channels and channel pairs use Julia's 1-based indices.

```julia
using DelayDifferentialAnalysis

result = run_dda_matrix(
    samples;
    device="cpu",
    channels=[1, 2, 3],
    flavors=["ST", "CT", "CD", "DE", "SY"],
    window_length=200,
    window_step=100,
    delays=[7, 10],
    model_terms=[1, 2, 10],
    derivative_points=4,
    order=4,
    nr_tau=2,
)

st = result["ST"]
println(st.matrix)
println(st.row_labels)
```

The CPU backend is the default. To use NVIDIA CUDA, install CUDA.jl in the
active environment and select a device with `device="cuda"` or
`device="cuda:0"`:

```julia
using Pkg
Pkg.add("CUDA")

result = run_dda_matrix(samples; device="cuda:0", flavors=["ST", "CT"])
```

CUDA.jl is optional and is loaded only for a CUDA device. Data preparation and
flavor assembly remain on the CPU; CUDA accelerates the independent regression
problems in batches. The native API and the external-binary API are separate,
so existing `run_DDA`, `run_st`, `run_ct`, and `run_de` calls are unchanged.

## Structure Selection

Structure selection has two recommended steps:

1. `structure_selection_compute` computes and caches DDA outputs.
2. `structure_selection_read` restores a previous compute run after restarting
   Julia.
3. `structure_selection_select` reads cached outputs and chooses the lowest-error
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
print_structure_selection(MOD)
```

The terminal printer shows `P_DDA`, a compact checkmark table for `MOD`, and a
Unicode model equation for each row. The polynomial order is inferred from the
number of `MOD` columns. `write_model_LaTeX` remains available when LaTeX output
is needed explicitly.

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
    num_cores=4,
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
for deterministic order. `num_cores` controls how many independent DDA binaries
run concurrently and defaults to `1`. Values above `Sys.CPU_THREADS` are capped
at the available logical CPU count.

The compute step writes `RUN1/structure_selection.toml` before starting the DDA
calls. Restore the same `StructureSelectionRun` in a later Julia session without
rerunning DDA:

```julia
run = structure_selection_read(prefix="RUN1")
```

### Select From Cached Outputs

```julia
selection = structure_selection_select(run; channels=[1])

println(selection.best_model)
println(selection.best_delays)
println(selection.best_score)
```

Use `models` to restrict selection to specific cached `MOD` row numbers:

```julia
selection = structure_selection_select(run; channels=[1], models=[1, 3, 7])
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

### Visualize Cached Selection Results

`plot_structure_selection` uses Plots.jl and reads the same cached ST outputs
without rerunning the DDA binary. The `:all` view shows the winning model number
for every valid two-delay pair after aggregating over all windows:

```julia
plot_structure_selection(run; mode=:all, channels=[1], models=[1, 3, 7])
```

The `:time` view selects the best model and delay pair independently at each
window and shows `tau1` and `tau2` against the raw DDA window coordinate `T`:

```julia
plot_structure_selection(run; mode=:time, channels=[1])
```

Both views compare only models with two active delays. Symmetric models are
mirrored across the delay-pair diagonal; asymmetric models retain their delay
ordering. Colors identify row numbers in `run.MOD`.

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
