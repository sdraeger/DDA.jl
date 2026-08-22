struct NativeFlavorResult
    id::String
    name::String
    matrix::Matrix{Float64}
    row_labels::Vector{String}
    window_markers::Vector{Float64}
end

struct NativeDDAResult
    flavors::Vector{NativeFlavorResult}
    window_markers::Vector{Float64}
    channel_labels::Vector{String}
end

function flavor_result(result::NativeDDAResult, flavor::AbstractString)::NativeFlavorResult
    id = _normalize_flavor(flavor)
    index = findfirst(item -> item.id == id, result.flavors)
    index === nothing && error("DDA flavor $id was not computed")
    return result.flavors[index]
end

Base.getindex(result::NativeDDAResult, flavor::AbstractString) = flavor_result(result, flavor)

struct ModelSpec
    derivative_points::Int
    window_length::Int
    window_step::Int
    max_delay::Int
    primary_terms::Vector{Vector{Int}}
end

struct PreparedWindow
    shifted::Matrix{Float64}
    derivative::Matrix{Float64}
    max_delay::Int
end

struct RegressionProblem
    design::Matrix{Float64}
    fit_target::Vector{Float64}
    residual_target::Vector{Float64}
end

struct SolvedBlock
    coefficients::Vector{Float64}
    rmse::Float64
end

_nan_block(feature_count::Int) = SolvedBlock(fill(NaN, feature_count), NaN)

struct WindowProblemRefs
    st::Vector{Int}
    ct::Vector{Int}
    de::Vector{Int}
    cd::Vector{Int}
    sy_forward::Vector{Int}
    sy_reverse::Vector{Int}
end
