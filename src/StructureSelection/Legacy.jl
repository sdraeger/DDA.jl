function _structure_selection(
    run_once::Function;
    file_path,
    channels=nothing,
    candidate_models=nothing,
    delays=nothing,
    candidate_delays=nothing,
    MOD=nothing,
    N_MOD=nothing,
    DDAorder=nothing,
    nr_delays::Integer=2,
    binary_path=nothing,
    derivative_points=nothing,
    order=nothing,
    tau_file=nothing,
    prefix=nothing,
    tau_file_suffix::AbstractString="",
    WL=nothing,
    WS=nothing,
    input_format=nothing,
    metric::Symbol=:mean_error,
    out_dir=nothing,
    model_scope=:joint,
    randomize::Bool=true,
    rng=Random.GLOBAL_RNG,
    cleanup_on_error::Bool=false,
    _trial_prefix::AbstractString="structure_selection",
    _artifacts_dir::Union{AbstractString, Nothing}=nothing,
    kwargs...,
)::Union{StructureSelectionResult, PerChannelStructureSelectionResult}
    haskey(kwargs, :tau_file_prefix) && error("Pass `prefix`, not `tau_file_prefix`")

    delay_spec = _resolve_structure_delays(delays, candidate_delays, tau_file)
    output_prefix = prefix === nothing ? nothing : expanduser(String(prefix))
    artifacts_dir = _resolve_structure_artifacts_dir(
        delay_spec.mode,
        out_dir,
        _artifacts_dir,
        output_prefix,
    )
    tau_prefix =
        delay_spec.mode == :explicit ? nothing : joinpath(String(artifacts_dir), "TAU_ALL__")
    created_artifacts_dir =
        delay_spec.mode == :pool && _artifacts_dir === nothing && output_prefix === nothing

    try
        return _structure_selection_resolved(
            run_once;
            file_path=file_path,
            channels=channels,
            candidate_models=candidate_models,
            delay_spec=delay_spec,
            MOD=MOD,
            N_MOD=N_MOD,
            DDAorder=DDAorder,
            nr_delays=nr_delays,
            binary_path=binary_path,
            derivative_points=derivative_points,
            order=order,
            tau_prefix=tau_prefix,
            tau_file_suffix=tau_file_suffix,
            WL=WL,
            WS=WS,
            input_format=input_format,
            metric=metric,
            out_dir=out_dir,
            model_scope=model_scope,
            randomize=randomize,
            rng=rng,
            cleanup_on_error=cleanup_on_error,
            _trial_prefix=_trial_prefix,
            _artifacts_dir=artifacts_dir,
            kwargs...,
        )
    catch
        if cleanup_on_error && created_artifacts_dir && artifacts_dir !== nothing
            rm(String(artifacts_dir); recursive=true, force=true)
        end
        rethrow()
    end
end

function _structure_selection_resolved(
    run_once::Function;
    file_path,
    channels,
    candidate_models,
    delay_spec,
    MOD,
    N_MOD,
    DDAorder,
    nr_delays::Integer,
    binary_path,
    derivative_points,
    order,
    tau_prefix,
    tau_file_suffix::AbstractString,
    WL,
    WS,
    input_format,
    metric::Symbol,
    out_dir,
    model_scope,
    randomize::Bool,
    rng,
    cleanup_on_error::Bool,
    _trial_prefix::AbstractString,
    _artifacts_dir::Union{AbstractString, Nothing},
    kwargs...,
)::Union{StructureSelectionResult, PerChannelStructureSelectionResult}
    scope = _normalize_model_scope(model_scope)
    if scope == :per_channel
        results = ChannelStructureSelectionResult[]
        for (channel_idx, channel) in enumerate(_normalize_structure_channels(channels))
            selection = _structure_selection_resolved(
                run_once;
                file_path=file_path,
                channels=[channel],
                candidate_models=candidate_models,
                delay_spec=delay_spec,
                MOD=MOD,
                N_MOD=N_MOD,
                DDAorder=DDAorder,
                nr_delays=nr_delays,
                binary_path=binary_path,
                derivative_points=derivative_points,
                order=order,
                tau_prefix=tau_prefix,
                tau_file_suffix=tau_file_suffix,
                WL=WL,
                WS=WS,
                input_format=input_format,
                metric=metric,
                out_dir=out_dir,
                model_scope=:joint,
                randomize=randomize,
                rng=rng,
                cleanup_on_error=cleanup_on_error,
                _trial_prefix="structure_selection_ch$(channel_idx)",
                _artifacts_dir=_artifacts_dir,
                kwargs...,
            )
            push!(results, ChannelStructureSelectionResult(channel_idx, channel, selection))
        end
        return PerChannelStructureSelectionResult(results)
    end

    derivative_points !== nothing || error("`derivative_points` is required for structure selection")
    model_order = _resolve_structure_order(order, DDAorder)

    models = _resolve_candidate_models(
        candidate_models;
        MOD=MOD,
        N_MOD=N_MOD,
        DDAorder=model_order,
        nr_delays=nr_delays,
    )
    P_DDA = _p_dda(model_order; nr_delays=nr_delays)
    output_root = delay_spec.mode == :pool ? String(_artifacts_dir) : _output_root(out_dir)

    trials = StructureSelectionTrial[]
    best_trial = nothing

    if delay_spec.mode == :explicit
        for model in models
            model_id = _model_filename_id(model, P_DDA; nr_delays=nr_delays, order=model_order)
            for (delay_idx, model_delays) in enumerate(delay_spec.values)
                out_fn = _trial_out_fn(output_root, model_id, delay_idx, _trial_prefix)
                result = run_once(;
                    file_path=file_path,
                    channels=channels,
                    flavors=["ST"],
                    binary_path=binary_path,
                    model=model,
                    delays=model_delays,
                    derivative_points=Int(derivative_points),
                    order=model_order,
                    nr_tau=length(model_delays),
                    WL=WL,
                    WS=WS,
                    input_format=input_format,
                    out_fn=out_fn,
                    kwargs...,
                )
                score = _score_result(result, metric)
                trial = StructureSelectionTrial(model, model_delays, score, result, out_fn)
                push!(trials, trial)
                if best_trial === nothing || trial.score < best_trial.score
                    best_trial = trial
                end
            end
        end
    else
        for model in _ordered_candidates(models, randomize, rng)
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
            result = _run_or_reuse_pool_output(out_fn, length(tau_rows)) do
                _write_tau_file(tau_path, tau_rows)
                run_once(;
                    file_path=file_path,
                    channels=channels,
                    flavors=["ST"],
                    binary_path=binary_path,
                    model=executable_model,
                    delays=first(tau_rows),
                    derivative_points=Int(derivative_points),
                    order=model_order,
                    nr_tau=nr,
                    tau_file=tau_path,
                    WL=WL,
                    WS=WS,
                    input_format=input_format,
                    out_fn=out_fn,
                    kwargs...,
                )
            end
            result === nothing && continue
            best_delays, score = _best_tau_row_score(result, metric, tau_rows)
            trial = StructureSelectionTrial(model, best_delays, score, result, out_fn, tau_path)
            push!(trials, trial)
            if best_trial === nothing || trial.score < best_trial.score
                best_trial = trial
            end
        end
    end

    best_trial === nothing && error("No usable structure-selection candidates were evaluated")
    return StructureSelectionResult(
        best_trial.model,
        best_trial.delays,
        best_trial.score,
        best_trial.result,
        trials,
        output_root,
    )
end
