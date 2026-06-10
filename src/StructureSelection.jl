"""Structure-selection utilities for DDA model and delay searches."""
module StructureSelection

using Statistics
using ..Runner: run_DDA

export StructureSelectionTrial, StructureSelectionResult, structure_selection

"""One evaluated structure-selection candidate."""
struct StructureSelectionTrial
    model::Any
    delays::Vector{Int}
    score::Float64
    result::Any
    out_fn::Union{String, Nothing}
end

"""Result returned by `structure_selection`."""
struct StructureSelectionResult
    best_model::Any
    best_delays::Vector{Int}
    best_score::Float64
    best_result::Any
    trials::Vector{StructureSelectionTrial}
end

"""
    structure_selection(; file_path, channels, binary_path=nothing,
        candidate_models, candidate_delays, derivative_points, order, kwargs...)

Evaluate each candidate model/delay combination with `run_DDA` and return the
candidate with the smallest ST error score. Candidate models use the same
encoding accepted by `run_DDA`: either `-MODEL` integer indices or model-matrix
rows that are converted to binary indices by the existing request path.
"""
function structure_selection(; kwargs...)::StructureSelectionResult
    return _structure_selection(run_DDA; kwargs...)
end

function _structure_selection(
    run_once::Function;
    file_path,
    channels,
    candidate_models,
    candidate_delays,
    binary_path=nothing,
    derivative_points=nothing,
    order=nothing,
    WL=nothing,
    WS=nothing,
    input_format=nothing,
    metric::Symbol=:mean_error,
    out_dir=nothing,
    kwargs...,
)::StructureSelectionResult
    derivative_points !== nothing || error("`derivative_points` is required for structure selection")
    order !== nothing || error("`order` is required for structure selection")

    models = _normalize_candidate_models(candidate_models)
    delay_sets = _normalize_candidate_delays(candidate_delays)
    output_root = _output_root(out_dir)

    trials = StructureSelectionTrial[]
    best_trial = nothing

    for (model_idx, model) in enumerate(models)
        for (delay_idx, delays) in enumerate(delay_sets)
            out_fn = _trial_out_fn(output_root, model_idx, delay_idx)
            result = run_once(;
                file_path=file_path,
                channels=channels,
                flavors=["ST"],
                binary_path=binary_path,
                model=model,
                delays=delays,
                derivative_points=Int(derivative_points),
                order=Int(order),
                nr_tau=length(delays),
                WL=WL,
                WS=WS,
                input_format=input_format,
                out_fn=out_fn,
                kwargs...,
            )
            score = _score_result(result, metric)
            trial = StructureSelectionTrial(model, delays, score, result, out_fn)
            push!(trials, trial)
            if best_trial === nothing || trial.score < best_trial.score
                best_trial = trial
            end
        end
    end

    best_trial === nothing && error("No structure-selection candidates were evaluated")
    return StructureSelectionResult(
        best_trial.model,
        best_trial.delays,
        best_trial.score,
        best_trial.result,
        trials,
    )
end

function _score_result(result, metric::Symbol)::Float64
    variant = _find_st_result(result)
    errors = Float64.(vec(getproperty(variant, :errors)))
    isempty(errors) && error("No ST error values found")

    if metric == :mean_error
        return mean(errors)
    elseif metric == :median_error
        return median(errors)
    elseif metric == :minimum_error
        return minimum(errors)
    end
    error("Unsupported structure-selection metric `$metric`")
end

function _find_st_result(result)
    for variant in getproperty(result, :variant_results)
        getproperty(variant, :variant_id) == "ST" && return variant
    end
    error("DDA result does not contain ST output")
end

function _normalize_candidate_models(candidate_models)::Vector{Any}
    if _is_model_candidate(candidate_models)
        return Any[candidate_models]
    end

    models = Any[]
    for model in candidate_models
        _is_model_candidate(model) || error(
            "`candidate_models` entries must be integer vectors or integer matrices",
        )
        push!(models, model)
    end
    isempty(models) && error("`candidate_models` must contain at least one model")
    return models
end

function _is_model_candidate(value)::Bool
    return value isa AbstractVector{<:Integer} || value isa AbstractMatrix{<:Integer}
end

function _normalize_candidate_delays(candidate_delays)::Vector{Vector{Int}}
    if candidate_delays isa AbstractVector{<:Integer}
        return [Int[candidate_delays...]]
    end

    delay_sets = Vector{Vector{Int}}()
    for delays in candidate_delays
        delays isa AbstractVector{<:Integer} || error(
            "`candidate_delays` entries must be integer vectors",
        )
        push!(delay_sets, Int[delays...])
    end
    isempty(delay_sets) && error("`candidate_delays` must contain at least one delay set")
    return delay_sets
end

function _output_root(out_dir)::Union{String, Nothing}
    out_dir === nothing && return nothing
    root = expanduser(String(out_dir))
    mkpath(root)
    return root
end

function _trial_out_fn(
    output_root::Union{String, Nothing},
    model_idx::Integer,
    delay_idx::Integer,
)::Union{String, Nothing}
    output_root === nothing && return nothing
    return joinpath(output_root, "structure_selection_m$(model_idx)_d$(delay_idx)")
end

end # module StructureSelection
