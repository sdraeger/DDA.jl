"""Structure-selection utilities for DDA model and delay searches."""
module StructureSelection

using Printf
using Random
using Statistics
using ..ModelEncoding: generate_monomials
using ..Runner: run_DDA

export ChannelStructureSelectionResult, PerChannelStructureSelectionResult
export StructureSelectionTrial, StructureSelectionResult, make_MOD, structure_selection
export StructureSelectionRun, structure_selection_compute, structure_selection_select
export print_structure_selection, write_model_terminal, write_model_LaTeX

const _POOL_LOCK_OWNER_FILE = "owner"
const _POOL_OWNERLESS_LOCK_GRACE_SECONDS = 60.0

"""One evaluated structure-selection candidate."""
struct StructureSelectionTrial
    model::Any
    delays::Vector{Int}
    score::Float64
    result::Any
    out_fn::Union{String, Nothing}
    tau_file::Union{String, Nothing}
end

StructureSelectionTrial(model, delays, score, result, out_fn) =
    StructureSelectionTrial(model, delays, score, result, out_fn, nothing)

"""Result returned by `structure_selection`."""
struct StructureSelectionResult
    best_model::Any
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
    channel_index::Int
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
explicit delay candidates.
Pass `prefix` to choose the output folder used for generated tau files and
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
        MOD=nothing, N_MOD=nothing, DDAorder, MOD_numbers=nothing, channels=nothing, ...)

Compute the cached DDA outputs needed for structure selection. This function
does not select a model.
"""
function structure_selection_compute(; kwargs...)::StructureSelectionRun
    return _structure_selection_compute(run_DDA; kwargs...)
end

"""
    structure_selection_select([run]; channels=nothing, channel=nothing,
        MOD_numbers=nothing, model_scope=:joint)

Select the best model and delays from cached structure-selection outputs.
This function only reads existing files.
"""
function structure_selection_select(run::StructureSelectionRun; kwargs...)
    return _structure_selection_select(run; kwargs...)
end

function structure_selection_select(; run=nothing, kwargs...)
    selected_run = run === nothing ? _LAST_STRUCTURE_SELECTION_RUN[] : run
    selected_run === nothing && error(
        "No structure-selection run is available; call `structure_selection_compute` or pass a `StructureSelectionRun`.",
    )
    return _structure_selection_select(selected_run; kwargs...)
end

"""
    make_MOD(N_MOD, DDAorder; nr_delays=2) -> Matrix{Int}

Generate the DDA structure-selection model library `MOD` used by the original
DDA scripts. `N_MOD` is the number, or list of numbers, of active monomial
terms per candidate model. `DDAorder` is the polynomial order used to construct
the `P_DDA` monomial table.
"""
function make_MOD(
    N_MOD,
    DDAorder::Integer;
    nr_delays::Integer=2,
)::Matrix{Int}
    return _structure_selection_model_space(N_MOD, DDAorder; nr_delays=nr_delays)
end

"""
    print_structure_selection([io], MOD, DDAorder; nr_delays=2, x="x")

Print the `P_DDA` monomial encoding table, `MOD` as a model-by-term checkmark
table, and each row of `MOD` as a Unicode model equation suitable for terminal
inspection.
"""
function print_structure_selection(
    io::IO,
    MOD::AbstractMatrix{<:Integer},
    DDAorder::Integer;
    nr_delays::Integer=2,
    x::AbstractString="x",
)::Nothing
    P_DDA = _p_dda(DDAorder; nr_delays=nr_delays)
    size(MOD, 2) == size(P_DDA, 1) || error(
        "MOD has $(size(MOD, 2)) columns, but P_DDA has $(size(P_DDA, 1)) monomials",
    )
    _validate_binary_MOD(MOD)

    println(io, "P_DDA")
    for idx in 1:size(P_DDA, 1)
        println(io, idx, ": ", collect(P_DDA[idx, :]))
    end

    println(io)
    _print_MOD_table(io, MOD)

    println(io)
    println(io, "Models")
    for mm in 1:size(MOD, 1)
        write_model_terminal(io, MOD, P_DDA, mm; x=x)
    end
    return nothing
end

function print_structure_selection(
    MOD::AbstractMatrix{<:Integer},
    DDAorder::Integer;
    kwargs...,
)::Nothing
    return print_structure_selection(stdout, MOD, DDAorder; kwargs...)
end

"""
    write_model_terminal([io], MOD, P_DDA, mm; x="x")

Print row `mm` of `MOD` as a Unicode equation for direct terminal inspection.
"""
function write_model_terminal(
    io::IO,
    MOD::AbstractMatrix{<:Integer},
    P_DDA::AbstractMatrix{<:Integer},
    mm::Integer;
    x::AbstractString="x",
)::Nothing
    1 <= mm <= size(MOD, 1) || error("Model row $mm out of range")
    model_terms = findall(==(1), vec(MOD[mm, :]))
    terms = P_DDA[model_terms, :]

    pieces = String[]
    for k in 1:size(terms, 1)
        monomial = _terminal_monomial(terms[k, :], x)
        push!(pieces, "a$(_subscript(k))·$monomial")
    end
    println(io, lpad(string(mm), 2), ": ẋ = ", join(pieces, " + "))
    return nothing
end

function write_model_terminal(
    MOD::AbstractMatrix{<:Integer},
    P_DDA::AbstractMatrix{<:Integer},
    mm::Integer;
    io::IO=stdout,
    kwargs...,
)::Nothing
    return write_model_terminal(io, MOD, P_DDA, mm; kwargs...)
end

"""
    write_model_LaTeX([io], MOD, SSYM, P_DDA, mm, x="x")

Print row `mm` of `MOD` as a LaTeX DDA equation using the supplied monomial
mapping table `P_DDA`.
"""
function write_model_LaTeX(
    io::IO,
    MOD::AbstractMatrix{<:Integer},
    _SSYM,
    P_DDA::AbstractMatrix{<:Integer},
    mm::Integer,
    x::AbstractString="x",
)::Nothing
    1 <= mm <= size(MOD, 1) || error("Model row $mm out of range")
    model_terms = findall(==(1), vec(MOD[mm, :]))
    terms = P_DDA[model_terms, :]
    nr_delays = maximum(P_DDA)

    @printf(io, "%2d & \\dot{%s} & = & ", mm, x)
    for k in 1:size(terms, 1)
        @printf(io, "a_%d ", k)
        for tau in 1:nr_delays
            exponent = count(==(tau), terms[k, :])
            if exponent == 1
                @printf(io, "%s_%d ", x, tau)
            elseif exponent > 1
                @printf(io, "%s_%d^%d ", x, tau, exponent)
            end
        end
        k < size(terms, 1) && print(io, "+ ")
    end
    println(io)
    return nothing
end

function write_model_LaTeX(
    MOD::AbstractMatrix{<:Integer},
    _SSYM,
    P_DDA::AbstractMatrix{<:Integer},
    mm::Integer,
    x::AbstractString="x";
    io::IO=stdout,
)::Nothing
    return write_model_LaTeX(io, MOD, _SSYM, P_DDA, mm, x)
end

function _terminal_monomial(row, x::AbstractString)::String
    parts = String[]
    nr_delays = maximum(row)
    for tau in 1:nr_delays
        exponent = count(==(tau), row)
        if exponent == 1
            push!(parts, "$x$(_subscript(tau))")
        elseif exponent > 1
            push!(parts, "$x$(_subscript(tau))$(_superscript(exponent))")
        end
    end
    return join(parts, "·")
end

function _subscript(value::Integer)::String
    return _translate_digits(value, ("₀", "₁", "₂", "₃", "₄", "₅", "₆", "₇", "₈", "₉"))
end

function _superscript(value::Integer)::String
    return _translate_digits(value, ("⁰", "¹", "²", "³", "⁴", "⁵", "⁶", "⁷", "⁸", "⁹"))
end

function _translate_digits(value::Integer, digits)::String
    value >= 0 || error("Unicode digit formatting only supports non-negative integers")
    return join(digits[Int(ch - '0') + 1] for ch in string(value))
end

function _print_MOD_table(io::IO, MOD::AbstractMatrix{<:Integer})::Nothing
    println(io, "MOD rows × P_DDA terms")
    println(io, "model | ", join(string.(1:size(MOD, 2)), " | "))
    for row_idx in 1:size(MOD, 1)
        cells = [MOD[row_idx, col_idx] == 1 ? "✓" : " " for col_idx in 1:size(MOD, 2)]
        println(io, row_idx, " | ", join(cells, " | "))
    end
    return nothing
end

function _validate_binary_MOD(MOD::AbstractMatrix{<:Integer})::Nothing
    all(value -> value == 0 || value == 1, MOD) || error("MOD must be a binary matrix")
    return nothing
end

function _structure_selection_model_space(
    N_MOD,
    DDAorder::Integer;
    nr_delays::Integer=2,
)::Matrix{Int}
    nr_delays == 2 || error("Structure-selection model generation currently supports nr_delays=2")
    DDAorder > 0 || error("DDAorder must be positive")

    P_DDA = _p_dda(DDAorder; nr_delays=nr_delays)
    model_sizes = _normalize_model_sizes(N_MOD)
    L = size(P_DDA, 1)
    mirror = _mirror_indices(P_DDA)

    rows = Vector{Vector{Int}}()
    for N in model_sizes
        1 <= N <= L || error("N_MOD entries must be between 1 and $L, got $N")
        combos = _combinations_indices(L, N)
        seen = Set{Tuple{Vararg{Int}}}()
        for combo in combos
            mirrored = Tuple(sort(mirror[combo]))
            key = min(Tuple(combo), mirrored)
            key in seen && continue
            push!(seen, key)

            row = zeros(Int, L)
            row[combo] .= 1
            push!(rows, row)
        end
    end

    isempty(rows) && error("No MOD rows generated")
    MOD = Matrix{Int}(undef, length(rows), L)
    for row_idx in eachindex(rows)
        MOD[row_idx, :] = rows[row_idx]
    end
    return MOD
end

function _p_dda(DDAorder::Integer; nr_delays::Integer=2)::Matrix{Int}
    nr_delays > 0 || error("nr_delays must be positive")
    DDAorder > 0 || error("DDAorder must be positive")
    monomials = generate_monomials(Int(nr_delays), Int(DDAorder))
    P_DDA = Matrix{Int}(undef, length(monomials), Int(DDAorder))
    for (row_idx, monomial) in enumerate(monomials)
        P_DDA[row_idx, :] = monomial
    end
    return P_DDA
end

function _normalize_model_sizes(N_MOD)::Vector{Int}
    if N_MOD isa Integer
        return [Int(N_MOD)]
    elseif N_MOD isa AbstractVector{<:Integer}
        sizes = Int.(collect(N_MOD))
        isempty(sizes) && error("N_MOD must contain at least one entry")
        return sizes
    end
    error("N_MOD must be an integer or vector/range of integers")
end

function _mirror_indices(P_DDA::AbstractMatrix{<:Integer})::Vector{Int}
    row_to_index = Dict(Tuple(P_DDA[idx, :]) => idx for idx in 1:size(P_DDA, 1))
    mirror = Vector{Int}(undef, size(P_DDA, 1))
    for idx in 1:size(P_DDA, 1)
        mirrored = _mirror_monomial(P_DDA[idx, :])
        mirror[idx] = row_to_index[Tuple(mirrored)]
    end
    return mirror
end

function _mirror_monomial(row)::Vector{Int}
    mirrored = [value == 1 ? 2 : value == 2 ? 1 : Int(value) for value in row]
    sort!(mirrored)
    return mirrored
end

function _combinations_indices(n::Integer, k::Integer)::Vector{Vector{Int}}
    result = Vector{Vector{Int}}()
    current = Int[]

    function visit(start::Int, remaining::Int)
        if remaining == 0
            push!(result, copy(current))
            return
        end
        last_start = Int(n) - remaining + 1
        for value in start:last_start
            push!(current, value)
            visit(value + 1, remaining - 1)
            pop!(current)
        end
    end

    visit(1, Int(k))
    return result
end

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
    _trial_prefix::AbstractString="structure_selection",
    kwargs...,
)::StructureSelectionRun
    haskey(kwargs, :tau_file_prefix) && error("Pass `prefix`, not `tau_file_prefix`")
    prefix === nothing && error("`prefix` is required for `structure_selection_compute`")
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

    for model_number in _ordered_model_numbers(model_numbers, randomize, rng)
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
        _run_or_reuse_pool_output(out_fn, length(tau_rows)) do
            _write_tau_file(tau_path, tau_rows)
            run_once(;
                file_path=file_path,
                channels=channel_list,
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
                load_results=false,
                kwargs...,
            )
        end
    end

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
    _LAST_STRUCTURE_SELECTION_RUN[] = run
    return run
end

function _structure_selection_select(
    run::StructureSelectionRun;
    channels=nothing,
    channel=nothing,
    MOD_numbers=nothing,
    model_scope=:joint,
    metric::Symbol=:mean_error,
)
    selected_channels = _resolve_selection_channel_argument(channels, channel)
    model_numbers = _resolve_selected_run_models(run, MOD_numbers)
    scope = _normalize_model_scope(model_scope)
    if scope == :per_channel
        groups = _selection_channel_groups(run, selected_channels, model_numbers, :per_channel)
        results = ChannelStructureSelectionResult[]
        for (idx, group) in enumerate(groups)
            selection = _structure_selection_select_joint(run, model_numbers, group, metric)
            push!(results, ChannelStructureSelectionResult(idx, first(group), selection))
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
        for model in _candidate_order(models, randomize, rng)
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

function _resolve_structure_delays(delays, candidate_delays, tau_file)
    tau_file === nothing || error(
        "`structure_selection` generates `-TAU_file` inputs from `delays`; pass `delays=...` and optionally `tau_file_suffix=...` instead of `tau_file`",
    )
    if delays !== nothing && candidate_delays !== nothing
        error("Pass `delays`, not both `delays` and deprecated `candidate_delays`")
    end
    selected = delays !== nothing ? delays : candidate_delays
    selected === nothing && error("Provide `delays` for structure selection")

    if selected isa AbstractVector{<:Integer}
        delay_pool = Int[selected...]
        isempty(delay_pool) && error("`delays` must contain at least one delay")
        return (mode=:pool, values=delay_pool)
    end

    return (mode=:explicit, values=_normalize_explicit_delay_sets(selected))
end

function _normalize_explicit_delay_sets(candidate_delays)::Vector{Vector{Int}}
    delay_sets = Vector{Vector{Int}}()
    for delays in candidate_delays
        delays isa AbstractVector{<:Integer} || error(
            "`delays` entries must be integer vectors when explicit delay candidates are supplied",
        )
        values = Int[delays...]
        isempty(values) && error("Explicit delay candidates cannot be empty")
        push!(delay_sets, values)
    end
    isempty(delay_sets) && error("`delays` must contain at least one delay candidate")
    return delay_sets
end

function _resolve_structure_artifacts_dir(mode::Symbol, out_dir, artifacts_dir, prefix)
    if mode == :explicit
        return artifacts_dir
    end
    artifacts_dir !== nothing && return String(artifacts_dir)
    if prefix !== nothing
        root = expanduser(String(prefix))
        mkpath(root)
        return root
    end

    parent = out_dir === nothing ? tempdir() : expanduser(String(out_dir))
    mkpath(parent)
    return mktempdir(parent; prefix="structure_selection_")
end

function _candidate_order(models::Vector{Any}, randomize::Bool, rng)::Vector{Any}
    randomize || return models
    length(models) <= 1 && return models
    return Random.shuffle(rng, models)
end

function _run_or_reuse_pool_output(
    run_candidate::Function,
    output_base::Union{String, Nothing},
    expected_tau_rows::Integer,
)
    output_base === nothing && return run_candidate()

    lock_path = "$(output_base).lock"
    while true
        existing = _existing_pool_result(output_base, expected_tau_rows)
        existing !== nothing && return existing

        if !ispath(lock_path) && _pool_artifact_conflict(output_base, expected_tau_rows)
            return nothing
        end

        if _try_create_lock(lock_path)
            try
                existing = _existing_pool_result(output_base, expected_tau_rows)
                existing !== nothing && return existing
                if _pool_artifact_conflict(output_base, expected_tau_rows)
                    return nothing
                end
                return run_candidate()
            finally
                rm(lock_path; recursive=true, force=true)
            end
        end

        existing = _wait_for_pool_result(output_base, lock_path, expected_tau_rows)
        existing !== nothing && return existing
    end
end

function _try_create_lock(lock_path::AbstractString)::Bool
    try
        _create_lock(lock_path)
        return true
    catch
        if ispath(lock_path)
            _remove_stale_lock(lock_path) || return false
            _create_lock(lock_path)
            return true
        end
        rethrow()
    end
end

function _create_lock(lock_path::AbstractString)::Nothing
    mkdir(lock_path)
    try
        _write_lock_owner(lock_path)
    catch
        rm(lock_path; recursive=true, force=true)
        rethrow()
    end
    return nothing
end

function _wait_for_pool_result(output_base::String, lock_path::String, expected_tau_rows::Integer)
    while ispath(lock_path)
        existing = _existing_pool_result(output_base, expected_tau_rows)
        existing !== nothing && return existing
        _remove_stale_lock(lock_path) && break
        sleep(0.5)
    end

    return _existing_pool_result(output_base, expected_tau_rows)
end

function _write_lock_owner(lock_path::AbstractString)::Nothing
    write(
        joinpath(lock_path, _POOL_LOCK_OWNER_FILE),
        "pid=$(getpid())\nhost=$(gethostname())\n",
    )
    return nothing
end

function _remove_stale_lock(lock_path::AbstractString)::Bool
    isdir(lock_path) || return false
    owner = _read_lock_owner(lock_path)
    if owner === nothing
        _ownerless_lock_is_stale(lock_path) || return false
        rm(lock_path; recursive=true, force=true)
        return true
    end
    if owner.host == gethostname() && !_pid_is_running(owner.pid)
        rm(lock_path; recursive=true, force=true)
        return true
    end
    return false
end

function _ownerless_lock_is_stale(lock_path::AbstractString, now::Real=time())::Bool
    return now - stat(lock_path).mtime > _POOL_OWNERLESS_LOCK_GRACE_SECONDS
end

function _read_lock_owner(lock_path::AbstractString)
    owner_path = joinpath(lock_path, _POOL_LOCK_OWNER_FILE)
    isfile(owner_path) || return nothing
    values = Dict{String, String}()
    try
        for line in eachline(owner_path)
            parts = split(line, "="; limit=2)
            length(parts) == 2 || continue
            values[String(parts[1])] = String(parts[2])
        end
        pid = parse(Int, get(values, "pid", "0"))
        host = get(values, "host", "")
        return (pid=pid, host=host)
    catch
        return nothing
    end
end

function _pid_is_running(pid::Integer)::Bool
    pid > 0 || return false
    Sys.iswindows() && return true
    result = ccall(:kill, Cint, (Cint, Cint), Cint(pid), Cint(0))
    return result == 0 || Base.Libc.errno() != 3
end

function _pool_artifact_conflict(output_base::String, expected_tau_rows::Integer)::Bool
    _existing_pool_result(output_base, expected_tau_rows) !== nothing && return false
    return isfile("$(output_base)_ST") || isfile("$(output_base).info")
end

function _existing_pool_result(output_base::String, expected_tau_rows::Integer)
    st_path = "$(output_base)_ST"
    isfile(st_path) || return nothing
    errors = _read_st_error_rows(st_path)
    errors === nothing && return nothing
    _has_complete_tau_rows(errors, expected_tau_rows) || return nothing
    return (variant_results=[(variant_id="ST", errors=errors)],)
end

function _has_complete_tau_rows(errors::AbstractMatrix, expected_tau_rows::Integer)::Bool
    expected_tau_rows > 0 || return false
    row_count = size(errors, 1)
    return row_count >= expected_tau_rows && row_count % expected_tau_rows == 0
end

function _read_st_error_rows(st_path::String)::Union{Matrix{Float64}, Nothing}
    values = Float64[]
    try
        for line in eachline(st_path)
            stripped = strip(line)
            isempty(stripped) && continue
            parts = split(stripped)
            isempty(parts) && continue
            push!(values, parse(Float64, parts[end]))
        end
    catch
        return nothing
    end
    isempty(values) && return nothing
    return reshape(values, :, 1)
end

function _normalize_model_scope(model_scope)::Symbol
    scope = Symbol(String(model_scope))
    (scope == :joint || scope == :per_channel) && return scope
    error("`model_scope` must be `:joint` or `:per_channel`, got `$model_scope`")
end

function _normalize_structure_channels(channels)::Vector{Int}
    channels === nothing && error("`channels` is required when `model_scope=:per_channel`")
    channels isa AbstractVector{<:Integer} || error(
        "`channels` must be an integer vector or range when `model_scope=:per_channel`",
    )
    channel_list = Int[channels...]
    isempty(channel_list) && error("`channels` must contain at least one channel")
    return channel_list
end

function _score_result(result, metric::Symbol)::Float64
    variant = _find_st_result(result)
    errors = Float64.(vec(getproperty(variant, :errors)))
    return _score_values(errors, metric)
end

function _find_st_result(result)
    for variant in getproperty(result, :variant_results)
        getproperty(variant, :variant_id) == "ST" && return variant
    end
    error("DDA result does not contain ST output")
end

function _resolve_structure_order(order, DDAorder)::Int
    order === nothing && DDAorder === nothing && error("`order` or `DDAorder` is required")
    if order !== nothing && DDAorder !== nothing && Int(order) != Int(DDAorder)
        error("`order` and `DDAorder` disagree")
    end
    selected = DDAorder === nothing ? order : DDAorder
    selected > 0 || error("DDAorder must be positive")
    return Int(selected)
end

function _resolve_candidate_models(
    candidate_models;
    MOD,
    N_MOD,
    DDAorder::Integer,
    nr_delays::Integer,
)::Vector{Any}
    if candidate_models !== nothing
        return _normalize_candidate_models(candidate_models)
    elseif MOD !== nothing
        return _models_from_MOD(MOD)
    elseif N_MOD !== nothing
        return _models_from_MOD(make_MOD(N_MOD, DDAorder; nr_delays=nr_delays))
    end
    error("Provide `candidate_models`, `MOD`, or `N_MOD` for structure selection")
end

function _resolve_MOD(MOD, N_MOD, DDAorder::Integer; nr_delays::Integer)::Matrix{Int}
    if MOD !== nothing
        MOD_matrix = Matrix{Int}(MOD)
        _validate_binary_MOD(MOD_matrix)
        P_DDA = _p_dda(DDAorder; nr_delays=nr_delays)
        size(MOD_matrix, 2) == size(P_DDA, 1) || error(
            "MOD has $(size(MOD_matrix, 2)) columns, but P_DDA has $(size(P_DDA, 1)) monomials",
        )
        return MOD_matrix
    end
    N_MOD !== nothing || error("Pass `MOD` or `N_MOD`")
    return make_MOD(N_MOD, DDAorder; nr_delays=nr_delays)
end

function _resolve_MOD_numbers(MOD::AbstractMatrix, MOD_numbers)::Vector{Int}
    numbers = MOD_numbers === nothing ? collect(1:size(MOD, 1)) : Int[MOD_numbers...]
    isempty(numbers) && error("`MOD_numbers` must contain at least one model row")
    all(number -> 1 <= number <= size(MOD, 1), numbers) ||
        error("`MOD_numbers` entries must be in 1:$(size(MOD, 1))")
    return numbers
end

function _resolve_selected_run_models(run::StructureSelectionRun, MOD_numbers)::Vector{Int}
    selected = MOD_numbers === nothing ? run.model_numbers : Int[MOD_numbers...]
    isempty(selected) && error("`MOD_numbers` must contain at least one model row")
    allowed = Set(run.model_numbers)
    all(number -> number in allowed, selected) || error(
        "`MOD_numbers` must refer to model rows computed in the supplied structure-selection run",
    )
    return selected
end

function _models_from_MOD(MOD::AbstractMatrix{<:Integer})::Vector{Any}
    _validate_binary_MOD(MOD)

    models = Any[]
    for row_idx in 1:size(MOD, 1)
        active = findall(==(1), vec(MOD[row_idx, :]))
        isempty(active) && error("MOD row $row_idx does not contain any active model terms")
        push!(models, active)
    end
    return models
end

function _model_from_MOD_row(MOD::AbstractMatrix{<:Integer}, model_number::Integer)::Vector{Int}
    1 <= model_number <= size(MOD, 1) || error("Model row $model_number out of range")
    active = findall(==(1), vec(MOD[model_number, :]))
    isempty(active) && error("MOD row $model_number does not contain any active model terms")
    return Int[active...]
end

function _ordered_model_numbers(numbers::Vector{Int}, randomize::Bool, rng)::Vector{Int}
    randomize || return numbers
    length(numbers) <= 1 && return numbers
    return Random.shuffle(rng, numbers)
end

function _model_symmetry(model, P_DDA::AbstractMatrix{<:Integer}; nr_delays::Integer, order::Integer)
    model_indices = _model_indices(model, P_DDA; nr_delays=nr_delays, order=order)
    terms = P_DDA[model_indices, :]
    positive_delays = unique(value for value in vec(terms) if value > 0)
    nr = length(positive_delays)
    1 <= nr <= 2 || error("Generated tau files currently support one or two active delay variables, got $nr")
    nr == 1 && return nr, 0

    mirrored_terms = Matrix{Int}(undef, size(terms)...)
    for row_idx in 1:size(terms, 1)
        mirrored_terms[row_idx, :] = _mirror_monomial(terms[row_idx, :])
    end
    sym = sortslices(Matrix{Int}(terms), dims=1) == sortslices(mirrored_terms, dims=1) ? 1 : 0
    return nr, sym
end

function _model_filename_id(model, P_DDA::AbstractMatrix{<:Integer}; nr_delays::Integer, order::Integer)::String
    indices = model isa AbstractVector{<:Integer} ? Int[model...] :
              _model_indices(model, P_DDA; nr_delays=nr_delays, order=order)
    return join(lpad.(string.(indices), 2, '0'), "_")
end

function _model_for_tau_file(
    model,
    P_DDA::AbstractMatrix{<:Integer};
    nr_delays::Integer,
    order::Integer,
    nr_tau::Integer,
)
    nr_tau == nr_delays && return model

    model_indices = _model_indices(model, P_DDA; nr_delays=nr_delays, order=order)
    terms = P_DDA[model_indices, :]
    active_delays = sort(unique(value for value in vec(terms) if value > 0))
    length(active_delays) == nr_tau || error(
        "Model uses $(length(active_delays)) active delay variables, but nr_tau=$nr_tau",
    )

    delay_map = Dict(delay => idx for (idx, delay) in enumerate(active_delays))
    remapped_terms = Matrix{Int}(undef, size(terms)...)
    for row_idx in 1:size(terms, 1)
        for col_idx in 1:size(terms, 2)
            value = Int(terms[row_idx, col_idx])
            remapped_terms[row_idx, col_idx] = value == 0 ? 0 : delay_map[value]
        end
    end

    executable_P_DDA = _p_dda(order; nr_delays=nr_tau)
    row_to_index = Dict(Tuple(executable_P_DDA[idx, :]) => idx for idx in 1:size(executable_P_DDA, 1))
    return Int[row_to_index[Tuple(remapped_terms[row_idx, :])] for row_idx in 1:size(remapped_terms, 1)]
end

function _model_indices(model, P_DDA::AbstractMatrix{<:Integer}; nr_delays::Integer, order::Integer)::Vector{Int}
    if model isa AbstractVector{<:Integer}
        indices = Int[model...]
        all(index -> 1 <= index <= size(P_DDA, 1), indices) || error("Model indices must be in 1:$(size(P_DDA, 1))")
        return indices
    elseif model isa AbstractMatrix{<:Integer}
        size(model, 2) == order || error("Model matrix must have $order columns")
        row_to_index = Dict(Tuple(P_DDA[idx, :]) => idx for idx in 1:size(P_DDA, 1))
        indices = Int[]
        for row_idx in 1:size(model, 1)
            row = Int[model[row_idx, col_idx] for col_idx in 1:size(model, 2)]
            any(value -> value < 0 || value > nr_delays, row) && error(
                "Model matrix row $row_idx has entries outside 0:$nr_delays: $(row)",
            )
            key = Tuple(row)
            haskey(row_to_index, key) || error("Model matrix row $row_idx does not match P_DDA: $(row)")
            push!(indices, row_to_index[key])
        end
        return indices
    end
    error("Model candidates must be integer vectors or matrices")
end

function _tau_rows(delay_pool::AbstractVector{<:Integer}, nr::Integer, sym::Integer)::Vector{Vector{Int}}
    delays = Int[delay_pool...]
    rows = Vector{Vector{Int}}()
    if nr == 1
        for tau in delays
            push!(rows, [tau])
        end
    elseif nr == 2 && sym == 0
        for tau1 in delays
            for tau2 in delays
                tau1 == tau2 && continue
                push!(rows, [tau1, tau2])
            end
        end
    elseif nr == 2 && sym == 1
        for idx1 in 1:length(delays)
            for idx2 in (idx1 + 1):length(delays)
                push!(rows, [delays[idx1], delays[idx2]])
            end
        end
    else
        error("Unsupported tau-file structure nr=$nr sym=$sym")
    end
    isempty(rows) && error("No tau rows generated for nr=$nr sym=$sym")
    return rows
end

function _tau_file_path(
    tau_prefix::AbstractString,
    nr::Integer,
    sym::Integer,
    suffix::AbstractString,
)::String
    return "$(tau_prefix)$(nr)_$(sym)$(suffix)"
end

function _read_tau_rows(tau_path::AbstractString)::Vector{Vector{Int}}
    rows = Vector{Vector{Int}}()
    for line in eachline(tau_path)
        stripped = strip(line)
        isempty(stripped) && continue
        push!(rows, parse.(Int, split(stripped)))
    end
    isempty(rows) && error("Tau file contains no rows: $tau_path")
    return rows
end

function _write_tau_file(
    tau_path::AbstractString,
    tau_rows::AbstractVector{<:AbstractVector{<:Integer}},
)::String
    artifacts_dir = dirname(tau_path)
    mkpath(artifacts_dir)
    contents = _tau_file_contents(tau_rows)
    isfile(tau_path) && read(tau_path, String) == contents && return String(tau_path)

    tmp_path, io = mktemp(artifacts_dir)
    try
        write(io, contents)
        close(io)
        mv(tmp_path, tau_path; force=true)
    catch
        isopen(io) && close(io)
        rm(tmp_path; force=true)
        rethrow()
    end
    return String(tau_path)
end

function _tau_file_contents(tau_rows::AbstractVector{<:AbstractVector{<:Integer}})::String
    buffer = IOBuffer()
    for row in tau_rows
        println(buffer, join(row, " "))
    end
    return String(take!(buffer))
end

function _resolve_selection_channel_argument(channels, channel)
    channels !== nothing && channel !== nothing && error("Pass `channels` or `channel`, not both")
    return channel === nothing ? channels : channel
end

function _selection_channel_groups(
    run::StructureSelectionRun,
    channels,
    model_numbers::AbstractVector{<:Integer},
    model_scope::Symbol,
)::Vector{Vector{Int}}
    if channels === nothing
        all_channels = run.channels === nothing ?
                       collect(1:_infer_cached_channel_count(run, model_numbers)) :
                       copy(run.channels)
        return model_scope == :per_channel ? [[channel] for channel in all_channels] : [all_channels]
    elseif channels isa AbstractVector{<:Integer}
        selected = Int[channels...]
        isempty(selected) && error("`channels` cannot be empty")
        return model_scope == :per_channel ? [[channel] for channel in selected] : [selected]
    end

    groups = Vector{Vector{Int}}()
    for group in channels
        group isa AbstractVector{<:Integer} || error("Channel groups must be integer vectors")
        values = Int[group...]
        isempty(values) && error("Channel groups cannot be empty")
        push!(groups, values)
    end
    isempty(groups) && error("`channels` cannot be empty")
    return groups
end

function _infer_cached_channel_count(
    run::StructureSelectionRun,
    model_numbers::AbstractVector{<:Integer},
)::Int
    P_DDA = _p_dda(run.DDAorder; nr_delays=run.nr_delays)
    for model_number in model_numbers
        model = _model_from_MOD_row(run.MOD, model_number)
        nr, sym = _model_symmetry(model, P_DDA; nr_delays=run.nr_delays, order=run.DDAorder)
        tau_path = _tau_file_path(joinpath(run.prefix, "TAU_ALL__"), nr, sym, run.tau_file_suffix)
        isfile(tau_path) || continue
        tau_rows = _read_tau_rows(tau_path)
        model_id = _model_filename_id(model, P_DDA; nr_delays=run.nr_delays, order=run.DDAorder)
        st_path = "$(_trial_out_fn(run.prefix, model_id, nothing, run.trial_prefix))_ST"
        errors = _read_structure_selection_errors(st_path, length(model) + 1, length(tau_rows))
        errors === nothing && continue
        return size(errors, 1)
    end
    error("Could not infer channel count from cached structure-selection outputs")
end

function _channel_positions(
    run::StructureSelectionRun,
    channels::AbstractVector{<:Integer},
    n_channels::Integer,
)::Vector{Int}
    if run.channels === nothing
        positions = Int[channels...]
    else
        positions = Int[]
        for channel in channels
            idx = findfirst(==(channel), run.channels)
            idx === nothing && error("Channel $channel was not part of the structure-selection compute run")
            push!(positions, idx)
        end
    end
    all(position -> 1 <= position <= n_channels, positions) ||
        error("Requested channels must resolve to positions in 1:$n_channels")
    return positions
end

function _read_structure_selection_errors(
    st_path::AbstractString,
    n_fields::Integer,
    n_tau_rows::Integer,
)::Union{Matrix{Float64}, Nothing}
    isfile(st_path) || return nothing
    raw = _read_numeric_matrix(st_path)
    raw === nothing && return nothing
    size(raw, 2) > 2 || return nothing
    payload = raw[:, 3:end]
    denominator = Int(n_fields) * Int(n_tau_rows)
    if denominator <= 0 || size(payload, 2) % denominator != 0
        return nothing
    end
    n_channels = div(size(payload, 2), denominator)
    n_channels > 0 || return nothing
    summary = Float64[median(view(payload, :, col)) for col in 1:size(payload, 2)]
    values = reshape(summary, Int(n_fields), n_channels, Int(n_tau_rows))
    return Matrix(values[end, :, :])
end

function _read_numeric_matrix(path::AbstractString)::Union{Matrix{Float64}, Nothing}
    rows = Vector{Vector{Float64}}()
    n_cols = 0
    try
        for line in eachline(path)
            stripped = strip(line)
            isempty(stripped) && continue
            values = parse.(Float64, split(stripped))
            if n_cols == 0
                n_cols = length(values)
            elseif length(values) != n_cols
                return nothing
            end
            push!(rows, values)
        end
    catch
        return nothing
    end
    isempty(rows) && return nothing
    matrix = Matrix{Float64}(undef, length(rows), n_cols)
    for row_idx in eachindex(rows)
        matrix[row_idx, :] = rows[row_idx]
    end
    return matrix
end

function _score_structure_error_rows(
    errors::AbstractMatrix{<:Real},
    channel_positions::AbstractVector{<:Integer},
    metric::Symbol,
)::Vector{Float64}
    return [
        _score_values(vec(errors[channel_positions, tau_idx]), metric)
        for tau_idx in 1:size(errors, 2)
    ]
end

function _best_tau_row_score(result, metric::Symbol, tau_rows::AbstractVector{<:AbstractVector{<:Integer}})
    variant = _find_st_result(result)
    errors = Float64.(getproperty(variant, :errors))
    scores = _score_error_rows(errors, metric, length(tau_rows))
    best_idx = argmin(scores)
    return Int[tau_rows[best_idx]...], scores[best_idx]
end

function _score_error_rows(errors, metric::Symbol, n_rows::Integer)::Vector{Float64}
    if errors isa AbstractVector && length(errors) == n_rows
        return Float64[errors...]
    elseif errors isa AbstractVector && length(errors) % n_rows == 0
        return [_score_values(errors[row_idx:n_rows:end], metric) for row_idx in 1:n_rows]
    elseif ndims(errors) == 2 && size(errors, 1) == n_rows
        return [_score_values(vec(errors[row_idx, :]), metric) for row_idx in 1:n_rows]
    elseif ndims(errors) == 2 && size(errors, 2) == n_rows
        return [_score_values(vec(errors[:, col_idx]), metric) for col_idx in 1:n_rows]
    elseif ndims(errors) == 2 && size(errors, 1) % n_rows == 0
        return [_score_values(vec(errors[row_idx:n_rows:end, :]), metric) for row_idx in 1:n_rows]
    end
    error("ST errors shape $(size(errors)) does not match generated tau-row count $n_rows")
end

function _score_values(values::AbstractVector{<:Real}, metric::Symbol)::Float64
    isempty(values) && error("No ST error values found")
    if metric == :mean_error
        return mean(Float64.(values))
    elseif metric == :median_error
        return median(Float64.(values))
    elseif metric == :minimum_error
        return minimum(Float64.(values))
    end
    error("Unsupported structure-selection metric `$metric`")
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

function _output_root(out_dir)::Union{String, Nothing}
    out_dir === nothing && return nothing
    root = expanduser(String(out_dir))
    mkpath(root)
    return root
end

function _trial_out_fn(
    output_root::Union{String, Nothing},
    model_id::AbstractString,
    delay_idx::Union{Integer, Nothing},
    prefix::AbstractString="structure_selection",
)::Union{String, Nothing}
    output_root === nothing && return nothing
    delay_suffix = delay_idx === nothing ? "" : "_d$(delay_idx)"
    return joinpath(output_root, "$(prefix)_$(model_id)$(delay_suffix)")
end

end # module StructureSelection
