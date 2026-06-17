"""Structure-selection utilities for DDA model and delay searches."""
module StructureSelection

using Printf
using Statistics
using ..ModelEncoding: generate_monomials
using ..Runner: run_DDA

export ChannelStructureSelectionResult, PerChannelStructureSelectionResult
export StructureSelectionTrial, StructureSelectionResult, make_MOD, structure_selection
export print_structure_selection, write_model_terminal, write_model_LaTeX

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

"""
    structure_selection(; file_path, channels, binary_path=nothing,
        candidate_models=nothing, MOD=nothing, N_MOD=nothing,
        candidate_delays=nothing, tau_file=nothing, derivative_points,
        order=nothing, DDAorder=nothing, model_scope=:joint, kwargs...)

Evaluate each candidate model/delay combination with `run_DDA` and return the
candidate with the smallest ST error score. Candidate models can be supplied
directly, as a `MOD` matrix from `make_MOD`, or by passing `N_MOD` plus
`DDAorder`. Delay candidates can be supplied directly or as rows in `tau_file`.
With `model_scope=:joint`, one model is selected across all channels. With
`model_scope=:per_channel`, one model is selected independently per channel.
"""
function structure_selection(; kwargs...)::Union{StructureSelectionResult, PerChannelStructureSelectionResult}
    return _structure_selection(run_DDA; kwargs...)
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
    channels,
    candidate_models=nothing,
    candidate_delays=nothing,
    MOD=nothing,
    N_MOD=nothing,
    DDAorder=nothing,
    nr_delays::Integer=2,
    binary_path=nothing,
    derivative_points=nothing,
    order=nothing,
    tau_file=nothing,
    WL=nothing,
    WS=nothing,
    input_format=nothing,
    metric::Symbol=:mean_error,
    out_dir=nothing,
    model_scope=:joint,
    _trial_prefix::AbstractString="structure_selection",
    kwargs...,
)::Union{StructureSelectionResult, PerChannelStructureSelectionResult}
    scope = _normalize_model_scope(model_scope)
    if scope == :per_channel
        results = ChannelStructureSelectionResult[]
        for (channel_idx, channel) in enumerate(_normalize_structure_channels(channels))
            selection = _structure_selection(
                run_once;
                file_path=file_path,
                channels=[channel],
                candidate_models=candidate_models,
                candidate_delays=candidate_delays,
                MOD=MOD,
                N_MOD=N_MOD,
                DDAorder=DDAorder,
                nr_delays=nr_delays,
                binary_path=binary_path,
                derivative_points=derivative_points,
                order=order,
                tau_file=tau_file,
                WL=WL,
                WS=WS,
                input_format=input_format,
                metric=metric,
                out_dir=out_dir,
                model_scope=:joint,
                _trial_prefix="structure_selection_ch$(channel_idx)",
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
    delay_sets = _resolve_candidate_delays(candidate_delays, tau_file)
    output_root = _output_root(out_dir)

    trials = StructureSelectionTrial[]
    best_trial = nothing

    for (model_idx, model) in enumerate(models)
        for (delay_idx, delays) in enumerate(delay_sets)
            out_fn = _trial_out_fn(output_root, model_idx, delay_idx, _trial_prefix)
            result = run_once(;
                file_path=file_path,
                channels=channels,
                flavors=["ST"],
                binary_path=binary_path,
                model=model,
                delays=delays,
                derivative_points=Int(derivative_points),
                order=model_order,
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

function _resolve_candidate_delays(candidate_delays, tau_file)::Vector{Vector{Int}}
    if candidate_delays !== nothing
        return _normalize_candidate_delays(candidate_delays)
    elseif tau_file !== nothing
        return _read_tau_file(tau_file)
    end
    error("Provide `candidate_delays` or `tau_file` for structure selection")
end

function _read_tau_file(path)::Vector{Vector{Int}}
    rows = Vector{Vector{Int}}()
    for (line_idx, line) in enumerate(readlines(expanduser(String(path))))
        stripped = strip(line)
        (isempty(stripped) || startswith(stripped, "#")) && continue
        values = Int[]
        for part in split(stripped)
            parsed = tryparse(Int, part)
            parsed !== nothing || error("Invalid integer in tau file line $line_idx: $part")
            push!(values, parsed)
        end
        isempty(values) || push!(rows, values)
    end
    isempty(rows) && error("Tau file contains no delay rows: $path")
    return rows
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
    prefix::AbstractString="structure_selection",
)::Union{String, Nothing}
    output_root === nothing && return nothing
    return joinpath(output_root, "$(prefix)_m$(model_idx)_d$(delay_idx)")
end

end # module StructureSelection
