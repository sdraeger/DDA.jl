function _resolve_structure_delays(delays, candidate_delays, tau_file)
    tau_file === nothing || error(
        "`structure_selection` generates `-TAU_file` inputs from `delays`; pass `delays=...` and optionally `tau_file_suffix=...` instead of `tau_file`",
    )
    if delays !== nothing && candidate_delays !== nothing
        error("Pass `delays`, not both `delays` and deprecated `candidate_delays`")
    end
    selected = delays !== nothing ? delays : candidate_delays
    if delays === nothing && candidate_delays !== nothing
        Base.depwarn(
            "`candidate_delays` is deprecated, use `delays`.",
            :_resolve_structure_delays,
        )
    end
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

function _ordered_candidates(values::Vector{T}, randomize::Bool, rng)::Vector{T} where {T}
    return randomize && length(values) > 1 ? Random.shuffle(rng, values) : values
end

# Locks guard against concurrent processes computing the same candidate into the
# same output base. A directory that can be created atomically acts as the lock;
# stale locks (crashed processes) are reclaimed after `_POOL_LOCK_STALE_SECONDS`.
# Duplicate compute after a steal is safe: outputs are content-checked before reuse.
const _POOL_LOCK_STALE_SECONDS = Ref(12 * 3600.0)

_pool_lock_path(output_base::AbstractString) = "$(output_base).lock"

function _lock_is_stale(lock_path::AbstractString, now::Real=time())::Bool
    return isdir(lock_path) && now - stat(lock_path).mtime > _POOL_LOCK_STALE_SECONDS[]
end

function _run_or_reuse_pool_output(
    run_candidate::Function,
    output_base::Union{String, Nothing},
    expected_tau_rows::Integer,
)
    output_base === nothing && return run_candidate()

    lock_path = _pool_lock_path(output_base)
    while true
        state, result = _pool_artifact_state(output_base, expected_tau_rows)
        state == :complete && return result
        state == :conflict && !ispath(lock_path) && return nothing

        if _acquire_pool_lock(lock_path)
            try
                state, result = _pool_artifact_state(output_base, expected_tau_rows)
                state == :complete && return result
                state == :conflict && return nothing
                return run_candidate()
            finally
                rm(lock_path; recursive=true, force=true)
            end
        end

        existing = _wait_for_pool_result(output_base, lock_path, expected_tau_rows)
        existing !== nothing && return existing
    end
end

function _acquire_pool_lock(lock_path::AbstractString)::Bool
    try
        mkdir(lock_path)
        return true
    catch
        if _lock_is_stale(lock_path)
            rm(lock_path; recursive=true, force=true)
            try
                mkdir(lock_path)
                return true
            catch
                return false
            end
        end
        return false
    end
end

function _wait_for_pool_result(output_base::String, lock_path::String, expected_tau_rows::Integer)
    while ispath(lock_path) && !_lock_is_stale(lock_path)
        state, result = _pool_artifact_state(output_base, expected_tau_rows)
        state == :complete && return result
        sleep(0.5)
    end

    state, result = _pool_artifact_state(output_base, expected_tau_rows)
    return state == :complete ? result : nothing
end

function _pool_artifact_state(output_base::String, expected_tau_rows::Integer)
    st_path = "$(output_base)_ST"
    if isfile(st_path)
        errors = _read_st_error_rows(st_path)
        if errors !== nothing && _has_complete_tau_rows(errors, expected_tau_rows)
            result = (variant_results=[(variant_id="ST", errors=errors)],)
            return :complete, result
        end
    end
    conflict = isfile(st_path) || isfile("$(output_base).info")
    return conflict ? (:conflict, nothing) : (:missing, nothing)
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
