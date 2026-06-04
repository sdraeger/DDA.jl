function run_st(;
    file_path::Union{AbstractString, Nothing}=nothing,
    data::Union{AbstractMatrix{<:Real}, Nothing}=nothing,
    channels::Union{AbstractVector{<:Integer}, Nothing}=nothing,
    kwargs...,
)::STResult
    if (file_path === nothing) == (data === nothing)
        error("Pass exactly one of `file_path` or `data`")
    end
    if file_path !== nothing
        channels !== nothing || error("`channels` is required when using `file_path`")
        return _run_st_file(file_path, channels; kwargs...)
    end
    channels === nothing || error("`channels` is not used with `data`; the matrix defines the channel set")
    return _run_st_matrix(data; kwargs...)
end

function run_ct(;
    file_path::Union{AbstractString, Nothing}=nothing,
    data::Union{AbstractMatrix{<:Real}, Nothing}=nothing,
    channels::Union{AbstractVector{<:Integer}, Nothing}=nothing,
    kwargs...,
)::CTResult
    if (file_path === nothing) == (data === nothing)
        error("Pass exactly one of `file_path` or `data`")
    end
    if file_path !== nothing
        channels !== nothing || error("`channels` is required when using `file_path`")
        return _run_ct_file(file_path, channels; kwargs...)
    end
    channels === nothing || error("`channels` is not used with `data`; the matrix defines the channel set")
    return _run_ct_matrix(data; kwargs...)
end

function run_de(;
    file_path::Union{AbstractString, Nothing}=nothing,
    data::Union{AbstractMatrix{<:Real}, Nothing}=nothing,
    channels::Union{AbstractVector{<:Integer}, Nothing}=nothing,
    kwargs...,
)::DEResult
    if (file_path === nothing) == (data === nothing)
        error("Pass exactly one of `file_path` or `data`")
    end
    if file_path !== nothing
        channels !== nothing || error("`channels` is required when using `file_path`")
        return _run_de_file(file_path, channels; kwargs...)
    end
    channels === nothing || error("`channels` is not used with `data`; the matrix defines the channel set")
    return _run_de_matrix(data; kwargs...)
end
