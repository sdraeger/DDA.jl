"""
    run_dda_matrix(samples; device="cpu", flavors=("ST",), ...)

Run the pure Julia DDA engine on a `samples × channels` matrix. Channel and
pair indices are 1-based. `device` accepts `"cpu"`, `"cuda"`, or `"cuda:N"`.
CUDA.jl is loaded only when a CUDA device is requested.
"""
function run_dda_matrix(
    samples::AbstractMatrix{<:Real};
    device::AbstractString="cpu",
    channels=nothing,
    flavors=("ST",),
    window_length::Int=NATIVE_WINDOW_LENGTH,
    window_step::Int=NATIVE_WINDOW_STEP,
    delays=NATIVE_DELAYS,
    model_terms=NATIVE_MODEL_TERMS,
    derivative_points::Int=NATIVE_DERIVATIVE_POINTS,
    order::Int=NATIVE_ORDER,
    nr_tau::Int=NATIVE_NR_TAU,
    ct_channel_pairs=nothing,
    cd_channel_pairs=nothing,
    ct_window_length=nothing,
    ct_window_step=nothing,
    channel_labels=nothing,
    start::Real=0.0,
    stop=nothing,
    normalization::AbstractString="zscore",
    nr_exclude::Int=10,
    derivative_step::Int=1,
)::NativeDDAResult
    data = _validate_samples(samples)
    row_count, channel_count = size(data)
    labels = _normalize_channel_labels(channel_labels, channel_count)
    selected_channels = _normalize_channels(channels, channel_count)
    enabled = Set(_normalize_flavor(String(flavor)) for flavor in flavors)
    isempty(enabled) && error("At least one DDA flavor must be enabled")
    unsupported = setdiff(enabled, Set(["ST", "CT", "CD", "DE", "SY"]))
    isempty(unsupported) || error("Unsupported Julia DDA flavor(s): $(join(sort(collect(unsupported)), ", "))")
    _parse_device(device)

    model = _build_model_spec(;
        window_length,
        window_step,
        delays,
        model_terms,
        derivative_points,
        order,
        nr_tau,
    )
    bounds_start, bounds_length = _analysis_bounds(start, stop, row_count)
    marker = model.window_length + model.max_delay + 2model.derivative_points
    required_rows = max(marker - 1, 0)
    bounds_length >= required_rows || error(
        "Selected range has $bounds_length samples but DDA requires at least $required_rows " *
        "samples (window_length + 2*derivative_points + max(delay) - 1)",
    )
    window_count = 1 + div(bounds_length - required_rows, model.window_step)
    markers = Float64[
        bounds_start + (window - 1) * model.window_step + marker
        for window in 1:window_count
    ]

    ct_groups = _resolve_ct_groups(
        selected_channels,
        ct_channel_pairs,
        ct_window_length,
        ct_window_step,
        channel_count,
    )
    de_groups = _sliding_groups(
        selected_channels,
        ct_window_length,
        ct_window_step,
    )
    cd_pairs = _resolve_cd_pairs(selected_channels, cd_channel_pairs, channel_count)
    sy_pairs = _sy_pairs(selected_channels)

    enabled_st = "ST" in enabled
    enabled_ct = "CT" in enabled
    enabled_cd = "CD" in enabled
    enabled_de = "DE" in enabled
    enabled_sy = "SY" in enabled && !isempty(sy_pairs)
    any((enabled_st, enabled_ct, enabled_cd, enabled_de, enabled_sy)) || error(
        "No DDA flavors are enabled for the selected channels",
    )

    st_matrix = enabled_st ? fill(NaN, length(selected_channels), window_count) : nothing
    ct_matrix = enabled_ct ? fill(NaN, length(ct_groups), window_count) : nothing
    cd_matrix = enabled_cd ? fill(NaN, length(cd_pairs), window_count) : nothing
    de_matrix = enabled_de ? fill(NaN, length(de_groups), window_count) : nothing
    sy_matrix = enabled_sy ? fill(NaN, length(sy_pairs), window_count) : nothing
    analysis_channels = _analysis_channels(selected_channels, ct_groups, de_groups, cd_pairs)

    jobs_per_window =
        ((enabled_st || enabled_cd || enabled_de) ? length(analysis_channels) : 0) +
        (enabled_ct ? length(ct_groups) : 0) +
        (enabled_de ? length(de_groups) : 0) +
        (enabled_cd ? length(cd_pairs) : 0) +
        (enabled_sy ? 2length(sy_pairs) : 0)
    windows_per_batch = max(1, min(32, div(2048, max(jobs_per_window, 1))))
    feature_count = length(model.primary_terms)

    for first_window in 1:windows_per_batch:window_count
        last_window = min(first_window + windows_per_batch - 1, window_count)
        problems = RegressionProblem[]
        references = WindowProblemRefs[]

        for window in first_window:last_window
            prepared = _prepare_window(
                data,
                bounds_start,
                model,
                window,
                String(normalization),
                nr_exclude,
                derivative_step,
            )
            st_refs = zeros(Int, channel_count)
            if enabled_st || enabled_cd || enabled_de
                for channel in analysis_channels
                    st_refs[channel] = _push_problem!(
                        problems,
                        _group_problem(prepared, [channel], model.primary_terms, model.window_length),
                    )
                end
            end
            ct_refs = enabled_ct ? Int[
                _push_problem!(
                    problems,
                    _group_problem(prepared, group, model.primary_terms, model.window_length),
                ) for group in ct_groups
            ] : Int[]
            de_refs = enabled_de ? Int[
                _push_problem!(
                    problems,
                    _group_problem(prepared, group, model.primary_terms, model.window_length),
                ) for group in de_groups
            ] : Int[]
            cd_refs = enabled_cd ? Int[
                _push_problem!(
                    problems,
                    _directed_problem(
                        prepared,
                        target,
                        source,
                        target,
                        model.primary_terms,
                        model.secondary_terms,
                        model.window_length,
                    ),
                ) for (target, source) in cd_pairs
            ] : Int[]
            sy_forward = enabled_sy ? Int[
                _push_problem!(
                    problems,
                    _directed_problem(
                        prepared,
                        left,
                        right,
                        right,
                        model.primary_terms,
                        model.secondary_terms,
                        model.window_length,
                    ),
                ) for (left, right) in sy_pairs
            ] : Int[]
            sy_reverse = enabled_sy ? Int[
                _push_problem!(
                    problems,
                    _directed_problem(
                        prepared,
                        right,
                        left,
                        left,
                        model.primary_terms,
                        model.secondary_terms,
                        model.window_length,
                    ),
                ) for (left, right) in sy_pairs
            ] : Int[]
            push!(references, WindowProblemRefs(st_refs, ct_refs, de_refs, cd_refs, sy_forward, sy_reverse))
        end

        solutions = _solve_problems(problems, device)
        for (offset, refs) in enumerate(references)
            window = first_window + offset - 1
            if st_matrix !== nothing
                for (row, channel) in enumerate(selected_channels)
                    block = _block(solutions, refs.st[channel], feature_count)
                    !isempty(block.coefficients) && (st_matrix[row, window] = block.coefficients[1])
                end
            end
            if ct_matrix !== nothing
                for row in eachindex(ct_groups)
                    block = _block(solutions, refs.ct[row], feature_count)
                    !isempty(block.coefficients) && (ct_matrix[row, window] = block.coefficients[1])
                end
            end
            if de_matrix !== nothing
                for row in eachindex(de_groups)
                    block = _block(solutions, refs.de[row], feature_count)
                    de_matrix[row, window] = _de_value(
                        de_groups[row], refs.st, solutions, block.rmse, feature_count
                    )
                end
            end
            if cd_matrix !== nothing
                for (row, (target, _)) in enumerate(cd_pairs)
                    baseline = _block(solutions, refs.st[target], feature_count).rmse
                    directed = _block(solutions, refs.cd[row], 2feature_count).rmse
                    cd_matrix[row, window] = _causal_improvement(baseline, directed)
                end
            end
            if sy_matrix !== nothing
                for row in eachindex(sy_pairs)
                    forward = _block(solutions, refs.sy_forward[row], 2feature_count).rmse
                    reverse = _block(solutions, refs.sy_reverse[row], 2feature_count).rmse
                    sy_matrix[row, window] = _synchronization_value(forward, reverse)
                end
            end
        end
    end

    results = NativeFlavorResult[]
    st_matrix !== nothing && push!(results, NativeFlavorResult(
        "ST", "Single Timeseries (ST)", st_matrix,
        _channel_labels(labels, selected_channels), markers,
    ))
    ct_matrix !== nothing && push!(results, NativeFlavorResult(
        "CT", "Cross-Timeseries (CT)", ct_matrix,
        _group_labels(labels, ct_groups, "&"), markers,
    ))
    cd_matrix !== nothing && push!(results, NativeFlavorResult(
        "CD", "Cross-Dynamical (CD)", cd_matrix,
        _pair_labels(labels, cd_pairs, " <- "), markers,
    ))
    de_matrix !== nothing && push!(results, NativeFlavorResult(
        "DE", "Dynamical Ergodicity (DE)", de_matrix,
        _group_labels(labels, de_groups, "&"), markers,
    ))
    sy_matrix !== nothing && push!(results, NativeFlavorResult(
        "SY", "Synchronization (SY)", sy_matrix,
        _pair_labels(labels, sy_pairs, " <-> "), markers,
    ))
    return NativeDDAResult(results, markers, labels)
end

function _sliding_groups(channels::Vector{Int}, window_length, window_step)::Vector{Vector{Int}}
    length_value = window_length === nothing ? length(channels) : Int(window_length)
    step = window_step === nothing ? max(length_value, 1) : max(Int(window_step), 1)
    (length_value <= 0 || length(channels) < length_value) && return Vector{Vector{Int}}()
    groups = Vector{Vector{Int}}()
    first = 1
    while first + length_value - 1 <= length(channels)
        push!(groups, copy(channels[first:(first + length_value - 1)]))
        first += step
    end
    return groups
end

function _normalize_groups(groups, channel_count::Int)::Vector{Vector{Int}}
    normalized = Vector{Vector{Int}}()
    for group in groups
        values = Int[group...]
        isempty(values) && error("Channel groups must not be empty")
        all(channel -> 1 <= channel <= channel_count, values) || error(
            "Channel group indices must be in 1:$channel_count",
        )
        push!(normalized, values)
    end
    return normalized
end

function _resolve_ct_groups(channels, pairs, window_length, window_step, channel_count)
    pairs !== nothing && !isempty(pairs) && return _normalize_groups(pairs, channel_count)
    return _sliding_groups(channels, window_length, window_step)
end

function _resolve_cd_pairs(channels, pairs, channel_count)::Vector{Tuple{Int,Int}}
    if pairs !== nothing && !isempty(pairs)
        groups = _normalize_groups(pairs, channel_count)
        all(group -> length(group) == 2, groups) || error("CD pairs must contain two channels")
        return Tuple{Int,Int}[(group[1], group[2]) for group in groups]
    end
    return Tuple{Int,Int}[
        (target, source) for target in channels for source in channels if target != source
    ]
end

function _sy_pairs(channels::Vector{Int})::Vector{Tuple{Int,Int}}
    return Tuple{Int,Int}[
        (channels[index], channels[index + 1]) for index in 1:2:(length(channels) - 1)
    ]
end

function _analysis_channels(st_channels, ct_groups, de_groups, cd_pairs)::Vector{Int}
    channels = Set(st_channels)
    foreach(group -> union!(channels, group), ct_groups)
    foreach(group -> union!(channels, group), de_groups)
    foreach(pair -> union!(channels, pair), cd_pairs)
    return sort(collect(channels))
end

_channel_labels(labels, channels) = String[labels[channel] for channel in channels]
_group_labels(labels, groups, joiner) = String[
    join(_channel_labels(labels, group), joiner) for group in groups
]
_pair_labels(labels, pairs, joiner) = String[
    "$(labels[left])$joiner$(labels[right])" for (left, right) in pairs
]
