function _run_from_file_or_data(
    file_runner::Function,
    ;
    file_path::Union{AbstractString, Nothing},
    data::Union{AbstractMatrix{<:Real}, Nothing},
    channels::Union{AbstractVector{<:Integer}, Nothing},
    kwargs...,
)
    if (file_path === nothing) == (data === nothing)
        error("Pass exactly one of `file_path` or `data`")
    end
    if file_path !== nothing
        channels !== nothing || error("`channels` is required when using `file_path`")
        return file_runner(file_path, channels; kwargs...)
    end
    channels === nothing || error("`channels` is not used with `data`; the matrix defines the channel set")
    return _run_matrix_with_temp(data, file_runner; kwargs...)
end

function run_st(;
    file_path::Union{AbstractString, Nothing}=nothing,
    data::Union{AbstractMatrix{<:Real}, Nothing}=nothing,
    channels::Union{AbstractVector{<:Integer}, Nothing}=nothing,
    kwargs...,
)::STResult
    return _run_from_file_or_data(_run_st_file;
        file_path=file_path, data=data, channels=channels, kwargs...)
end

function run_ct(;
    file_path::Union{AbstractString, Nothing}=nothing,
    data::Union{AbstractMatrix{<:Real}, Nothing}=nothing,
    channels::Union{AbstractVector{<:Integer}, Nothing}=nothing,
    kwargs...,
)::CTResult
    return _run_from_file_or_data(_run_ct_file;
        file_path=file_path, data=data, channels=channels, kwargs...)
end

function run_de(;
    file_path::Union{AbstractString, Nothing}=nothing,
    data::Union{AbstractMatrix{<:Real}, Nothing}=nothing,
    channels::Union{AbstractVector{<:Integer}, Nothing}=nothing,
    kwargs...,
)::DEResult
    return _run_from_file_or_data(_run_de_file;
        file_path=file_path, data=data, channels=channels, kwargs...)
end
