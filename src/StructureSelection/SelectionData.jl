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
    isempty(selected) && error("`models`/`MOD_numbers` must contain at least one model row")
    allowed = Set(run.model_numbers)
    all(number -> number in allowed, selected) || error(
        "`models`/`MOD_numbers` must refer to model rows computed in the supplied structure-selection run",
    )
    return selected
end

function _models_from_MOD(MOD::AbstractMatrix{<:Integer})::Vector{Any}
    _validate_binary_MOD(MOD)
    return Any[_model_from_MOD_row(MOD, row_idx) for row_idx in axes(MOD, 1)]
end

function _model_from_MOD_row(MOD::AbstractMatrix{<:Integer}, model_number::Integer)::Vector{Int}
    1 <= model_number <= size(MOD, 1) || error("Model row $model_number out of range")
    active = findall(==(1), vec(MOD[model_number, :]))
    isempty(active) && error("MOD row $model_number does not contain any active model terms")
    return Int[active...]
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
        return model_matrix_to_encoding(
            model;
            num_delays=Int(nr_delays),
            polynomial_order=Int(order),
        )
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
