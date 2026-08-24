function _structure_selection_plot_data(
    run::StructureSelectionRun;
    mode::Symbol=:all,
    channels=nothing,
    models=nothing,
    metric::Symbol=:mean_error,
)
    mode in (:all, :time) || error("`mode` must be `:all` or `:time`")
    model_numbers = _resolve_selected_run_models(run, models)
    groups = _selection_channel_groups(run, channels, model_numbers, :joint)
    candidates = _two_delay_plot_candidates(run, model_numbers, first(groups))
    return mode == :all ?
           _structure_selection_all_data(run, candidates, metric) :
           _structure_selection_time_data(candidates, metric)
end

function _two_delay_plot_candidates(
    run::StructureSelectionRun,
    model_numbers::AbstractVector{<:Integer},
    channels::AbstractVector{<:Integer},
)
    P_DDA = _p_dda(run.DDAorder; nr_delays=run.nr_delays)
    candidates = []
    for model_number in model_numbers
        model = _model_from_MOD_row(run.MOD, model_number)
        candidate = _pool_candidate(
            model,
            P_DDA,
            run.delays;
            nr_delays=run.nr_delays,
            order=run.DDAorder,
            tau_prefix=joinpath(run.prefix, "TAU_ALL__"),
            output_root=run.prefix,
            tau_file_suffix=run.tau_file_suffix,
            trial_prefix=run.trial_prefix,
        )
        candidate.nr == 2 || continue
        isfile(candidate.tau_path) || continue
        tau_rows = _read_tau_rows(candidate.tau_path)
        all(row -> length(row) == 2, tau_rows) || error("Expected delay pairs in $(candidate.tau_path)")
        series = _read_structure_selection_series("$(candidate.out_fn)_ST", length(model) + 1, length(tau_rows))
        series === nothing && continue
        positions = _channel_positions(run, channels, size(series.errors, 1))
        push!(candidates, (
            model_number=Int(model_number),
            symmetry=candidate.sym,
            tau_rows=tau_rows,
            T=vec(series.T[:, 1]),
            errors=series.errors[positions, :, :],
        ))
    end
    isempty(candidates) && error("No usable cached outputs for two-delay models were found")
    return candidates
end

function _structure_selection_all_data(run::StructureSelectionRun, candidates, metric::Symbol)
    delays = sort(unique(run.delays))
    delay_indices = Dict(delay => idx for (idx, delay) in enumerate(delays))
    model_numbers = zeros(Int, length(delays), length(delays))
    best_scores = fill(Inf, length(delays), length(delays))

    for candidate in candidates
        scores = _aggregate_candidate_errors(candidate.errors, metric)
        for (tau_idx, pair) in enumerate(candidate.tau_rows)
            _set_delay_pair_winner!(
                model_numbers,
                best_scores,
                delay_indices,
                pair[1],
                pair[2],
                candidate.model_number,
                scores[tau_idx],
            )
            if candidate.symmetry == 1
                _set_delay_pair_winner!(
                    model_numbers,
                    best_scores,
                    delay_indices,
                    pair[2],
                    pair[1],
                    candidate.model_number,
                    scores[tau_idx],
                )
            end
        end
    end
    any(>(0), model_numbers) || error("No finite structure-selection scores were found")
    return (
        mode=:all,
        tau1=delays,
        tau2=delays,
        model_numbers=model_numbers,
    )
end

function _set_delay_pair_winner!(
    model_numbers::Matrix{Int},
    best_scores::Matrix{Float64},
    delay_indices::Dict{Int,Int},
    tau1::Integer,
    tau2::Integer,
    model_number::Integer,
    score::Real,
)
    isfinite(score) || return
    haskey(delay_indices, tau1) && haskey(delay_indices, tau2) || error(
        "Cached delay pair ($tau1, $tau2) is outside the run's delay pool",
    )
    row = delay_indices[Int(tau2)]
    col = delay_indices[Int(tau1)]
    current_model = model_numbers[row, col]
    if score < best_scores[row, col] ||
       (score == best_scores[row, col] && (current_model == 0 || model_number < current_model))
        model_numbers[row, col] = Int(model_number)
        best_scores[row, col] = Float64(score)
    end
    return
end

function _aggregate_candidate_errors(errors::AbstractArray{<:Real,3}, metric::Symbol)
    n_channels, n_tau, _ = size(errors)
    scores = Vector{Float64}(undef, n_tau)
    channel_scores = Vector{Float64}(undef, n_channels)
    for tau_idx in 1:n_tau
        for channel_idx in 1:n_channels
            channel_scores[channel_idx] = median(view(errors, channel_idx, tau_idx, :))
        end
        scores[tau_idx] = _score_values(channel_scores, metric)
    end
    return scores
end

function _structure_selection_time_data(candidates, metric::Symbol)
    T = copy(first(candidates).T)
    all(candidate -> candidate.T == T, candidates) || error(
        "Cached structure-selection outputs do not share the same window coordinates",
    )
    model_numbers = zeros(Int, length(T))
    tau1 = zeros(Int, length(T))
    tau2 = zeros(Int, length(T))
    best_scores = fill(Inf, length(T))

    for candidate in candidates
        scores = _time_candidate_errors(candidate.errors, metric)
        for window_idx in eachindex(T)
            for (tau_idx, pair) in enumerate(candidate.tau_rows)
                score = scores[tau_idx, window_idx]
                _is_better_time_winner(
                    score,
                    candidate.model_number,
                    pair,
                    best_scores[window_idx],
                    model_numbers[window_idx],
                    tau1[window_idx],
                    tau2[window_idx],
                ) || continue
                best_scores[window_idx] = score
                model_numbers[window_idx] = candidate.model_number
                tau1[window_idx], tau2[window_idx] = pair
            end
        end
    end
    any(>(0), model_numbers) || error("No finite structure-selection scores were found")
    return (
        mode=:time,
        T=T,
        tau1=tau1,
        tau2=tau2,
        model_numbers=model_numbers,
    )
end

function _time_candidate_errors(errors::AbstractArray{<:Real,3}, metric::Symbol)
    _, n_tau, n_windows = size(errors)
    scores = Matrix{Float64}(undef, n_tau, n_windows)
    for window_idx in 1:n_windows
        for tau_idx in 1:n_tau
            scores[tau_idx, window_idx] = _score_values(
                view(errors, :, tau_idx, window_idx),
                metric,
            )
        end
    end
    return scores
end

function _is_better_time_winner(
    score::Real,
    model_number::Integer,
    pair::AbstractVector{<:Integer},
    best_score::Real,
    best_model::Integer,
    best_tau1::Integer,
    best_tau2::Integer,
)::Bool
    isfinite(score) || return false
    score < best_score && return true
    score == best_score || return false
    best_model == 0 && return true
    return isless(
        (Int(model_number), Int(pair[1]), Int(pair[2])),
        (Int(best_model), Int(best_tau1), Int(best_tau2)),
    )
end
