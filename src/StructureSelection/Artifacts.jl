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

function _ordered_candidates(values::Vector{T}, randomize::Bool, rng)::Vector{T} where {T}
    return randomize && length(values) > 1 ? Random.shuffle(rng, values) : values
end

function _run_or_reuse_pool_output(
    run_candidate::Function,
    output_base::Union{String, Nothing},
    expected_tau_rows::Integer,
)
    output_base === nothing && return run_candidate()

    lock_path = "$(output_base).lock"
    while true
        state, result = _pool_artifact_state(output_base, expected_tau_rows)
        state == :complete && return result
        state == :conflict && !ispath(lock_path) && return nothing

        if _try_create_lock(lock_path)
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
        state, result = _pool_artifact_state(output_base, expected_tau_rows)
        state == :complete && return result
        _remove_stale_lock(lock_path) && break
        sleep(0.5)
    end

    state, result = _pool_artifact_state(output_base, expected_tau_rows)
    return state == :complete ? result : nothing
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
