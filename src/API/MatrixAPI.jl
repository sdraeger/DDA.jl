"""Run a file-based API function on an in-memory channels × samples matrix."""
function _run_matrix_with_temp(
    data::AbstractMatrix{<:Real},
    file_runner::Function;
    kwargs...,
)
    path = _write_temp_ascii(data)
    try
        return file_runner(path, collect(axes(data, 1)); kwargs...)
    finally
        isfile(path) && rm(path; force=true)
    end
end
