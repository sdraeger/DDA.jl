# =============================================================================
# STRUCTURED OUTPUT TYPES
# =============================================================================

"""A single timepoint's parsed data for one channel.

For stride-1 flavors (DE/SY) the single output value is stored in `value`
(it is an ergodicity/synchronization measure, not an error).
"""
struct StructuredTimepoint
    window_start::Int64
    window_end::Int64
    coefficients::Vector{Float64}
    value::Float64
end

"""All timepoints for a single channel/pair."""
struct StructuredChannelData
    channel_index::Int
    timepoints::Vector{StructuredTimepoint}
end

# =============================================================================
# LEGACY RESULT TYPES (backward compat)
# =============================================================================

"""
Result data for a single variant.

`T` contains the first two integer columns emitted by the binary (this differs
from `STResult.T`, which holds only the first column). `A` contains every
remaining column emitted by the binary, preserving output-file order.
"""
struct VariantResultData
    variant_id::String
    variant_name::String
    A::Matrix{Float64}
    coefficients::Array{Float64,3}
    errors::Matrix{Float64}
    T::Matrix{Int64}
    t::Vector{Float64}
    window_starts::Vector{Int64}
    window_ends::Vector{Int64}
    channel_labels::Union{Vector{String}, Nothing}
end

"""
DDA analysis result.

Flavors are accessed via [`flavor_matrix`](@ref) or indexing
(`result["ST"]`, `result["CT"]`, ...), each returning that flavor's raw output
matrix. Dynamic `result.ST`-style property access still works but is deprecated.
Top-level `T`, `t`, and `A` mirror the primary variant for backward compatibility.
"""
struct DDAResult
    id::String
    file_path::String
    channels::Vector{String}
    T::Matrix{Int64}
    t::Vector{Float64}
    A::Matrix{Float64}
    variant_results::Vector{VariantResultData}
    window_params::WindowParameters
    delay_params::DelayParameters
    created_at::Dates.DateTime
end

function DDAResult(
    id::String,
    file_path::String,
    channels::Vector{String},
    T::Matrix{Int64},
    A::Matrix{Float64},
    variant_results::Vector{VariantResultData},
    window_params::WindowParameters,
    delay_params::DelayParameters,
    created_at::Dates.DateTime,
)
    t = isempty(variant_results) ? Float64[] : first(variant_results).t
    return DDAResult(
        id,
        file_path,
        channels,
        T,
        t,
        A,
        variant_results,
        window_params,
        delay_params,
        created_at,
    )
end

const DDA_RESULT_FLAVOR_PROPERTIES = (:ST, :CT, :CD, :DE, :SY)

function _variant_result_property(result::DDAResult, name::Symbol)::Union{VariantResultData, Nothing}
    for variant in getfield(result, :variant_results)
        Symbol(getfield(variant, :variant_id)) == name && return variant
    end
    return nothing
end

function Base.getproperty(result::DDAResult, name::Symbol)
    name in fieldnames(typeof(result)) && return getfield(result, name)
    if name in DDA_RESULT_FLAVOR_PROPERTIES
        variant = _variant_result_property(result, name)
        if variant !== nothing
            Base.depwarn(
                "`result.$name` is deprecated; use `flavor_matrix(result, \"$name\")`.",
                :getproperty,
            )
            return getfield(variant, :A)
        end
        error("DDAResult does not contain flavor `$name`")
    end
    return getfield(result, name)
end

"""
    flavor_matrix(result::DDAResult, flavor) -> Matrix{Float64}

Return the raw output matrix of a flavor (`"ST"`, `"CT"`, `"CD"`, `"DE"`, `"SY"`).
This replaces the deprecated dynamic `result.ST`-style property access.
"""
function flavor_matrix(result::DDAResult, flavor::AbstractString)::Matrix{Float64}
    variant = _variant_result_property(result, Symbol(flavor))
    variant === nothing && error("DDAResult does not contain flavor `$flavor`")
    return getfield(variant, :A)
end

Base.getindex(result::DDAResult, flavor::AbstractString) = flavor_matrix(result, flavor)

function Base.propertynames(result::DDAResult, private::Bool=false)
    base_names = collect(fieldnames(typeof(result)))
    for variant in getfield(result, :variant_results)
        name = Symbol(getfield(variant, :variant_id))
        if !(name in base_names)
            push!(base_names, name)
        end
    end
    return Tuple(base_names)
end
