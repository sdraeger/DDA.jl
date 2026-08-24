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

"""
    run_st(; file_path, channels, data=nothing, kwargs...) -> STResult

Run single-timeseries (ST) DDA on an EDF/ASCII file using 1-based channel
indices, or directly on a `data` matrix (channels × samples).

Important keywords:
- `binary_path`: resolve the DDA binary without environment variables
- `WL`, `WS`: analysis window length and step (`-WL`, `-WS`)
- `delays`: delay (tau) values, default `[7, 10]`
- `model`: optional custom model; a vector of `-MODEL` indices or a monomial
  matrix (passing a custom `model` requires explicit `derivative_points`
  and `order`)
- `derivative_points`: `-dm` parameter (the `dm` alias is deprecated)
- `sampling_rate`: scalar or two-tuple passed to `-SR`; also scales `result.t`
- `TM`: offset used only for `result.t`; defaults to `max(delays)`
- `channel_labels`: optional labels; inferred from EDF/ASCII headers otherwise
- `out_fn`: keep raw output files at this base path instead of a temporary one

Returns an [`STResult`](@ref) with `coefficients`, `errors`, `T`, `t`,
`window_starts`/`window_ends`, `channel_labels`, and `params`.
`to_dataframe(result)` converts it to a long-format table.
"""
function run_st(;
    file_path::Union{AbstractString, Nothing}=nothing,
    data::Union{AbstractMatrix{<:Real}, Nothing}=nothing,
    channels::Union{AbstractVector{<:Integer}, Nothing}=nothing,
    kwargs...,
)::STResult
    return _run_from_file_or_data(_run_st_file;
        file_path=file_path, data=data, channels=channels, kwargs...)
end

"""
    run_ct(; file_path, channels, data=nothing, kwargs...) -> CTResult

Run cross-timeseries (CT) DDA on an EDF/ASCII file using 1-based channel
indices, or directly on a `data` matrix (channels × samples). Every channel
pair `(i, j)` with `i < j` is analyzed; the binary runs once per pair.

Extra keywords beyond `run_st`:
- `ct_wl`, `ct_ws`: CT-specific window length/step forwarded as `-WL_CT`/
  `-WS_CT`; also used for the result's window bounds

Requires at least two channels. Returns a [`CTResult`](@ref) with one row per
pair and `pair_labels` formatted as `"a-b"`.
"""
function run_ct(;
    file_path::Union{AbstractString, Nothing}=nothing,
    data::Union{AbstractMatrix{<:Real}, Nothing}=nothing,
    channels::Union{AbstractVector{<:Integer}, Nothing}=nothing,
    kwargs...,
)::CTResult
    return _run_from_file_or_data(_run_ct_file;
        file_path=file_path, data=data, channels=channels, kwargs...)
end

"""
    run_de(; file_path, channels, data=nothing, kwargs...) -> DEResult

Run dynamical-ergodicity (DE) DDA on an EDF/ASCII file using 1-based channel
indices, or directly on a `data` matrix (channels × samples). DE requires the
CT window parameters; its result window bounds derive from `ct_wl`/`ct_ws`.

Returns a [`DEResult`](@ref) with one `ergodicity` value per window.
"""
function run_de(;
    file_path::Union{AbstractString, Nothing}=nothing,
    data::Union{AbstractMatrix{<:Real}, Nothing}=nothing,
    channels::Union{AbstractVector{<:Integer}, Nothing}=nothing,
    kwargs...,
)::DEResult
    return _run_from_file_or_data(_run_de_file;
        file_path=file_path, data=data, channels=channels, kwargs...)
end
