# =============================================================================
# STRUCTURED OUTPUT TYPES
# =============================================================================

"""A single timepoint's parsed data for one channel."""
struct StructuredTimepoint
    window_start::Float64
    window_end::Float64
    coefficients::Vector{Float64}
    error::Float64
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

`T` contains the first two integer columns emitted by the binary. `A` contains every
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

Returned flavors are exposed as dynamic properties, for example `result.ST`
and `result.CT`, with each flavor property returning that flavor's `A` matrix
directly. Top-level `T`, `t`, and `A` mirror the primary variant for backward
compatibility.
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
    created_at::String
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
    created_at::String,
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
        variant !== nothing && return getfield(variant, :A)
        error("DDAResult does not contain flavor `$name`")
    end
    return getfield(result, name)
end

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
