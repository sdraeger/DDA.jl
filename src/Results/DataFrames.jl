# DataFrame conversion with a lazy DataFrames.jl import.
# invokelatest below is required: DataFrames methods are defined after
# this module precompiles (world-age).

import ..OptionalDeps

_dataframes_module()::Module = OptionalDeps.require(:DataFrames)

"""
    to_dataframe(r::STResult)

Convert an STResult to a long-format DataFrame.

Columns: `channel`, `window_start`, `window_end`, `a_1`, ..., `a_N`, `error`.

Requires DataFrames.jl to be installed.
"""
function to_dataframe(r::STResult)
    _dataframes_module()
    Base.invokelatest(_to_dataframe_st, r)
end

function _to_dataframe_st(r::STResult)
    return _to_dataframe_coeff_result(r, :channel, r.channel_labels)
end

function _to_dataframe_coeff_result(
    r::Union{STResult,CTResult},
    entity_column::Symbol,
    labels::Vector{String},
)
    DF = _dataframes_module()
    n_entities = size(r.coefficients, 1)
    nw = n_windows(r)
    nc = n_coeffs(r)

    rows = n_entities * nw
    col_entity = Vector{String}(undef, rows)
    col_T = Vector{Float64}(undef, rows)
    col_t = Vector{Float64}(undef, rows)
    col_wstart = Vector{Int64}(undef, rows)
    col_wend = Vector{Int64}(undef, rows)
    coeff_cols = [Vector{Float64}(undef, rows) for _ in 1:nc]
    col_error = Vector{Float64}(undef, rows)

    idx = 0
    for entity in 1:n_entities, w in 1:nw
        idx += 1
        col_entity[idx] = labels[entity]
        col_T[idx] = r.T[w]
        col_t[idx] = r.t[w]
        col_wstart[idx] = r.window_starts[w]
        col_wend[idx] = r.window_ends[w]
        for c in 1:nc
            coeff_cols[c][idx] = r.coefficients[entity, w, c]
        end
        col_error[idx] = r.errors[entity, w]
    end

    pairs = [entity_column => col_entity,
             Symbol("T") => col_T,
             Symbol("t") => col_t,
             Symbol("window_start") => col_wstart,
             Symbol("window_end") => col_wend]
    for c in 1:nc
        push!(pairs, Symbol("a_$c") => coeff_cols[c])
    end
    push!(pairs, Symbol("error") => col_error)

    return Base.invokelatest(DF.DataFrame, pairs...)
end

"""
    to_dataframe(r::CTResult)

Convert a CTResult to a long-format DataFrame.
"""
function to_dataframe(r::CTResult)
    _dataframes_module()
    Base.invokelatest(_to_dataframe_ct, r)
end

function _to_dataframe_ct(r::CTResult)
    return _to_dataframe_coeff_result(r, :pair, r.pair_labels)
end

"""
    to_dataframe(r::DEResult)

Convert a DEResult to a long-format DataFrame.
"""
function to_dataframe(r::DEResult)
    _dataframes_module()
    Base.invokelatest(_to_dataframe_de, r)
end

function _to_dataframe_de(r::DEResult)
    DF = _dataframes_module()
    nw = n_windows(r)

    col_T = r.T
    col_t = r.t
    col_wstart = r.window_starts
    col_wend = r.window_ends
    col_erg = r.ergodicity

    return Base.invokelatest(
        DF.DataFrame,
        Symbol("T") => col_T,
        Symbol("t") => col_t,
        Symbol("window_start") => col_wstart,
        Symbol("window_end") => col_wend,
        Symbol("ergodicity") => col_erg,
    )
end
