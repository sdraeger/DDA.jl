function _structure_selection_compute(
    run_once::Function;
    file_path,
    prefix,
    channels=nothing,
    MOD=nothing,
    N_MOD=nothing,
    DDAorder=nothing,
    order=nothing,
    nr_delays::Integer=2,
    delays=nothing,
    candidate_delays=nothing,
    derivative_points,
    binary_path=nothing,
    input_format=nothing,
    WL=nothing,
    WS=nothing,
    MOD_numbers=nothing,
    tau_file_suffix::AbstractString="",
    randomize::Bool=true,
    rng=Random.GLOBAL_RNG,
    num_cores::Integer=1,
    _trial_prefix::AbstractString="structure_selection",
    kwargs...,
)::StructureSelectionRun
    haskey(kwargs, :tau_file_prefix) && error("Pass `prefix`, not `tau_file_prefix`")
    prefix === nothing && error("`prefix` is required for `structure_selection_compute`")
    num_cores >= 1 || error("`num_cores` must be at least 1")
    delay_spec = _resolve_structure_delays(delays, candidate_delays, nothing)
    delay_spec.mode == :pool || error("`structure_selection_compute` expects `delays` as a flat delay pool")
    derivative_points !== nothing || error("`derivative_points` is required")
    model_order = _resolve_structure_order(order, DDAorder)
    MOD_matrix = _resolve_MOD(MOD, N_MOD, model_order; nr_delays=nr_delays)
    model_numbers = _resolve_MOD_numbers(MOD_matrix, MOD_numbers)
    output_root = expanduser(String(prefix))
    mkpath(output_root)
    P_DDA = _p_dda(model_order; nr_delays=nr_delays)
    tau_prefix = joinpath(output_root, "TAU_ALL__")
    channel_list = channels === nothing ? nothing : Int[channels...]
    run = StructureSelectionRun(
        output_root,
        MOD_matrix,
        model_order,
        Int(nr_delays),
        Int[delay_spec.values...],
        channel_list,
        model_numbers,
        Int(derivative_points),
        String(tau_file_suffix),
        String(_trial_prefix),
    )
    _write_structure_selection_run(run)

    candidates = map(_ordered_candidates(model_numbers, randomize, rng)) do model_number
        model = _model_from_MOD_row(MOD_matrix, model_number)
        nr, sym = _model_symmetry(model, P_DDA; nr_delays=nr_delays, order=model_order)
        tau_rows = _tau_rows(delay_spec.values, nr, sym)
        tau_path = _tau_file_path(tau_prefix, nr, sym, tau_file_suffix)
        executable_model = _model_for_tau_file(
            model,
            P_DDA;
            nr_delays=nr_delays,
            order=model_order,
            nr_tau=nr,
        )
        model_id = _model_filename_id(model, P_DDA; nr_delays=nr_delays, order=model_order)
        out_fn = _trial_out_fn(output_root, model_id, nothing, _trial_prefix)
        return (;
            model=executable_model,
            nr,
            tau_rows,
            tau_path,
            out_fn,
        )
    end

    for candidate in candidates
        _write_tau_file(candidate.tau_path, candidate.tau_rows)
    end

    run_candidate = function(candidate)
        _run_or_reuse_pool_output(candidate.out_fn, length(candidate.tau_rows)) do
            run_once(;
                file_path=file_path,
                channels=channel_list,
                flavors=["ST"],
                binary_path=binary_path,
                model=candidate.model,
                delays=first(candidate.tau_rows),
                derivative_points=Int(derivative_points),
                order=model_order,
                nr_tau=candidate.nr,
                tau_file=candidate.tau_path,
                WL=WL,
                WS=WS,
                input_format=input_format,
                out_fn=candidate.out_fn,
                load_results=false,
                kwargs...,
            )
        end
    end

    worker_count = _structure_selection_worker_count(num_cores, length(candidates))
    if worker_count == 1
        foreach(run_candidate, candidates)
    else
        asyncmap(run_candidate, candidates; ntasks=worker_count)
    end

    _LAST_STRUCTURE_SELECTION_RUN[] = run
    return run
end

function _structure_selection_worker_count(num_cores::Integer, candidate_count::Integer)::Int
    num_cores >= 1 || error("`num_cores` must be at least 1")
    candidate_count >= 1 || error("Structure selection requires at least one model")
    return min(Int(num_cores), Sys.CPU_THREADS, Int(candidate_count))
end

function _structure_selection_select(
    run::StructureSelectionRun;
    channels=nothing,
    channel=nothing,
    models=nothing,
    MOD_numbers=nothing,
    model_scope=:joint,
    metric::Symbol=:mean_error,
)
    selected_channels = _resolve_selection_channel_argument(channels, channel)
    models !== nothing && MOD_numbers !== nothing && error(
        "Pass `models` or `MOD_numbers`, not both",
    )
    selected_models = models === nothing ? MOD_numbers : models
    model_numbers = _resolve_selected_run_models(run, selected_models)
    scope = _normalize_model_scope(model_scope)
    if scope == :per_channel
        groups = _selection_channel_groups(run, selected_channels, model_numbers, :per_channel)
        results = ChannelStructureSelectionResult[]
        for group in groups
            selection = _structure_selection_select_joint(run, model_numbers, group, metric)
            push!(results, ChannelStructureSelectionResult(first(group), selection))
        end
        return PerChannelStructureSelectionResult(results)
    end

    groups = _selection_channel_groups(run, selected_channels, model_numbers, :joint)
    return _structure_selection_select_joint(run, model_numbers, first(groups), metric)
end

function _structure_selection_select_joint(
    run::StructureSelectionRun,
    model_numbers::AbstractVector{<:Integer},
    channels::AbstractVector{<:Integer},
    metric::Symbol,
)::StructureSelectionResult
    P_DDA = _p_dda(run.DDAorder; nr_delays=run.nr_delays)
    trials = StructureSelectionTrial[]
    best_trial = nothing

    for model_number in model_numbers
        model = _model_from_MOD_row(run.MOD, model_number)
        nr, sym = _model_symmetry(model, P_DDA; nr_delays=run.nr_delays, order=run.DDAorder)
        tau_rows = _tau_rows(run.delays, nr, sym)
        tau_path = _tau_file_path(joinpath(run.prefix, "TAU_ALL__"), nr, sym, run.tau_file_suffix)
        isfile(tau_path) || continue
        tau_rows = _read_tau_rows(tau_path)
        model_id = _model_filename_id(model, P_DDA; nr_delays=run.nr_delays, order=run.DDAorder)
        out_fn = _trial_out_fn(run.prefix, model_id, nothing, run.trial_prefix)
        st_path = "$(out_fn)_ST"
        errors = _read_structure_selection_errors(st_path, length(model) + 1, length(tau_rows))
        errors === nothing && continue
        positions = _channel_positions(run, channels, size(errors, 1))
        scores = _score_structure_error_rows(errors, positions, metric)
        best_idx = argmin(scores)
        result = (variant_results=[(variant_id="ST", errors=errors)],)
        trial = StructureSelectionTrial(model, Int[tau_rows[best_idx]...], scores[best_idx], result, out_fn, tau_path)
        push!(trials, trial)
        if best_trial === nothing || trial.score < best_trial.score
            best_trial = trial
        end
    end

    best_trial === nothing && error("No usable cached structure-selection outputs were found")
    return StructureSelectionResult(
        best_trial.model,
        best_trial.delays,
        best_trial.score,
        best_trial.result,
        trials,
        run.prefix,
    )
end
