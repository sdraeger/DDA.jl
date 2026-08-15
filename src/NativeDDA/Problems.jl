function _evaluate_term(
    prepared::PreparedWindow,
    channel::Int,
    sample::Int,
    delays::Vector{Int},
)::Float64
    product = 1.0
    for delay in delays
        row = sample + max(prepared.max_delay - delay, 0)
        value = prepared.shifted[row, channel]
        isnan(value) && return NaN
        product *= value
    end
    return product
end

function _group_problem(
    prepared::PreparedWindow,
    channels::Vector{Int},
    terms::Vector{Vector{Int}},
    window_length::Int,
)::Union{RegressionProblem,Nothing}
    total_rows = length(channels) * window_length
    total_rows == 0 && return nothing
    design = Matrix{Float64}(undef, total_rows, length(terms))
    target = Vector{Float64}(undef, total_rows)
    valid_rows = 0

    for channel in channels, sample in 1:window_length
        target_value = prepared.derivative[channel, sample]
        isnan(target_value) && continue
        values = Float64[_evaluate_term(prepared, channel, sample, term) for term in terms]
        any(isnan, values) && continue
        valid_rows += 1
        design[valid_rows, :] = values
        target[valid_rows] = target_value
    end
    valid_rows / total_rows >= 0.60 || return nothing
    values = copy(target[1:valid_rows])
    return RegressionProblem(copy(design[1:valid_rows, :]), values, copy(values))
end

function _directed_problem(
    prepared::PreparedWindow,
    primary_channel::Int,
    secondary_channel::Int,
    response_channel::Int,
    primary_terms::Vector{Vector{Int}},
    secondary_terms::Vector{Vector{Int}},
    window_length::Int,
)::Union{RegressionProblem,Nothing}
    feature_count = length(primary_terms) + length(secondary_terms)
    design = Matrix{Float64}(undef, window_length, feature_count)
    fit_target = Vector{Float64}(undef, window_length)
    residual_target = Vector{Float64}(undef, window_length)
    valid_rows = 0

    for sample in 1:window_length
        fit_value = prepared.derivative[primary_channel, sample]
        isnan(fit_value) && continue
        secondary = Float64[
            _evaluate_term(prepared, secondary_channel, sample, term)
            for term in secondary_terms
        ]
        any(isnan, secondary) && continue
        primary = Float64[
            _evaluate_term(prepared, primary_channel, sample, term)
            for term in primary_terms
        ]
        any(isnan, primary) && continue
        valid_rows += 1
        design[valid_rows, :] = vcat(secondary, primary)
        fit_target[valid_rows] = fit_value
        residual_target[valid_rows] = prepared.derivative[response_channel, sample]
    end
    valid_rows / window_length >= 0.60 || return nothing
    return RegressionProblem(
        copy(design[1:valid_rows, :]),
        copy(fit_target[1:valid_rows]),
        copy(residual_target[1:valid_rows]),
    )
end

function _push_problem!(
    problems::Vector{RegressionProblem},
    problem::Union{RegressionProblem,Nothing},
)::Int
    problem === nothing && return 0
    push!(problems, problem)
    return length(problems)
end

function _block(solutions::Vector{SolvedBlock}, reference::Int, feature_count::Int)::SolvedBlock
    return iszero(reference) ? _nan_block(feature_count) : solutions[reference]
end

function _de_value(
    channels::Vector{Int},
    st_references::Vector{Int},
    solutions::Vector{SolvedBlock},
    ct_rmse::Float64,
    feature_count::Int,
)::Float64
    (isempty(channels) || isnan(ct_rmse) || iszero(ct_rmse)) && return NaN
    baseline = mean(
        _block(solutions, st_references[channel], feature_count).rmse
        for channel in channels
    )
    return isnan(baseline) ? NaN : abs(baseline / ct_rmse - 1.0)
end

function _causal_improvement(baseline::Float64, causal::Float64)::Float64
    return isnan(baseline) || isnan(causal) ? NaN : baseline - causal
end

function _synchronization_value(forward::Float64, reverse::Float64)::Float64
    (isnan(forward) || isnan(reverse) || iszero(forward) || iszero(reverse)) && return NaN
    return reverse / forward - forward / reverse
end
