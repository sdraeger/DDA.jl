"""Structure-selection utilities for DDA model and delay searches."""
module StructureSelection

using Printf
using Random
using Statistics
using TOML
using ..ModelEncoding: generate_monomials, model_matrix_to_encoding
using ..Runner
using ..Runner: run_DDA

"""A model candidate: `-MODEL` indices, or a monomial matrix row-wise encoded."""
const ModelCandidate = Union{AbstractVector{Int}, AbstractMatrix{Int}}

"""One evaluated structure-selection candidate."""
struct StructureSelectionTrial
    model::ModelCandidate
    delays::Vector{Int}
    score::Float64
    # Live DDAResults from the one-shot API, cached-output stubs from select.
    result::Any
    out_fn::Union{String, Nothing}
    tau_file::Union{String, Nothing}
end

StructureSelectionTrial(model, delays, score, result, out_fn) =
    StructureSelectionTrial(model, delays, score, result, out_fn, nothing)

"""Result returned by `structure_selection`."""
struct StructureSelectionResult
    best_model::ModelCandidate
    best_delays::Vector{Int}
    best_score::Float64
    best_result::Any
    trials::Vector{StructureSelectionTrial}
    artifacts_dir::Union{String, Nothing}
end

StructureSelectionResult(best_model, best_delays, best_score, best_result, trials) =
    StructureSelectionResult(best_model, best_delays, best_score, best_result, trials, nothing)

"""Structure-selection result for one input channel in per-channel mode."""
struct ChannelStructureSelectionResult
    channel::Int
    selection::StructureSelectionResult
end

"""Result returned by `structure_selection(...; model_scope=:per_channel)`."""
struct PerChannelStructureSelectionResult
    results::Vector{ChannelStructureSelectionResult}
end

"""Metadata for a structure-selection DDA output cache."""
struct StructureSelectionRun
    prefix::String
    MOD::Matrix{Int}
    DDAorder::Int
    nr_delays::Int
    delays::Vector{Int}
    channels::Union{Vector{Int}, Nothing}
    model_numbers::Vector{Int}
    derivative_points::Int
    tau_file_suffix::String
    trial_prefix::String
end

const _LAST_STRUCTURE_SELECTION_RUN = Ref{Union{StructureSelectionRun, Nothing}}(nothing)

"""
    structure_selection(; file_path, channels=nothing, binary_path=nothing,
        candidate_models=nothing, MOD=nothing, N_MOD=nothing,
        delays=nothing, candidate_delays=nothing, derivative_points,
        order=nothing, DDAorder=nothing, model_scope=:joint, prefix=nothing,
        randomize=true,
        kwargs...)

Evaluate each candidate model/delay combination with `run_DDA` and return the
candidate with the smallest ST error score. Candidate models can be supplied
directly, as a `MOD` matrix from `make_MOD`, or by passing `N_MOD` plus
`DDAorder`. Pass `delays` as a flat delay pool, such as `(derivative_points + 1):TM`,
to generate Claudia-style `TAU_ALL__...` files, or as nested vectors for
explicit delay candidates. `candidate_delays` is a deprecated alias for
`delays`. Pass `prefix` to choose the output folder used for generated tau files and
structure-selection outputs, for example `/scratch/run42`.
Pool-mode candidates are evaluated in randomized order by default so concurrent
structure-selection runs are less likely to start with the same model.
With `model_scope=:joint`, one model is selected across all channels. With
`model_scope=:per_channel`, one model is selected independently per channel.
"""
function structure_selection(; kwargs...)::Union{StructureSelectionResult, PerChannelStructureSelectionResult}
    return _structure_selection(run_DDA; kwargs...)
end

"""
    structure_selection_compute(; file_path, prefix, delays, derivative_points,
        MOD=nothing, N_MOD=nothing, DDAorder, MOD_numbers=nothing,
        channels=nothing, num_cores=1, ...)

Compute the cached DDA outputs needed for structure selection. This function
does not select a model. Independent model calls run concurrently when
`num_cores > 1`, capped at `Sys.CPU_THREADS`.
"""
function structure_selection_compute(; kwargs...)::StructureSelectionRun
    return _structure_selection_compute(run_DDA; kwargs...)
end

"""
    structure_selection_read(prefix)
    structure_selection_read(; prefix)

Read a `StructureSelectionRun` from a result folder created by
`structure_selection_compute`. No DDA processes are started.
"""
function structure_selection_read(prefix::AbstractString)::StructureSelectionRun
    return _structure_selection_read(prefix)
end

structure_selection_read(; prefix)::StructureSelectionRun = structure_selection_read(prefix)

"""
    structure_selection_select([run]; channels=nothing, channel=nothing,
        models=nothing, MOD_numbers=nothing, model_scope=:joint)

Select the best model and delays from cached structure-selection outputs.
`models` contains the cached `MOD` row numbers to compare. If omitted, every
model computed in `run` is considered. This function only reads existing files.
"""
function structure_selection_select(run::StructureSelectionRun; kwargs...)
    return _structure_selection_select(run; kwargs...)
end

function structure_selection_select(; run=nothing, kwargs...)
    if run === nothing
        Base.depwarn(
            "Calling `structure_selection_select()` without `run` is deprecated; " *
            "pass the `StructureSelectionRun` returned by `structure_selection_compute` " *
            "or `structure_selection_read`.",
            :structure_selection_select,
        )
        selected_run = _LAST_STRUCTURE_SELECTION_RUN[]
        selected_run === nothing && error(
            "No structure-selection run is available; call `structure_selection_compute` or pass a `StructureSelectionRun`.",
        )
        return _structure_selection_select(selected_run; kwargs...)
    end
    return _structure_selection_select(run; kwargs...)
end


include("StructureSelection/Models.jl")
include("StructureSelection/Artifacts.jl")
include("StructureSelection/SelectionData.jl")
include("StructureSelection/Visualization.jl")
include("StructureSelection/RunMetadata.jl")
include("StructureSelection/Compute.jl")
include("StructureSelection/Legacy.jl")

end # module StructureSelection
