using Dates

struct Bounds
    start::Int
    stop::Int
end

struct AlgorithmSelection
    enabled_variants::Vector{String}
    select_mask::Union{String,Nothing}
end

struct WindowParameters
    window_length::Int
    window_step::Int
    ct_window_length::Union{Int,Nothing}
    ct_window_step::Union{Int,Nothing}
end

struct DelayParameters
    delays::Vector{Int}
end

struct DDARequest
    file_path::String
    channels::Union{Vector{Int},Nothing}
    bounds::Union{Bounds,Nothing}
    algorithm_selection::AlgorithmSelection
    window_parameters::WindowParameters
    delay_parameters::DelayParameters
    ct_channel_pairs::Union{Vector{Vector{Int}},Nothing}
end

struct VariantResult
    variant_id::String
    variant_name::String
    q_matrix::Matrix{Float64}
    channel_labels::Union{Vector{String},Nothing}
end

struct DDAResult
    file_path::String
    channels::Vector{String}
    q_matrix::Matrix{Float64}
    variant_results::Union{Vector{VariantResult},Nothing}
    raw_output::Union{String,Nothing}
    window_parameters::WindowParameters
    delay_parameters::DelayParameters
end

function DDAResult(
    file_path::String,
    channels::Vector{String},
    q_matrix::Matrix{Float64},
    window_parameters::WindowParameters,
    delay_parameters::DelayParameters
)
    DDAResult(
        file_path,
        channels,
        q_matrix,
        nothing,
        nothing,
        window_parameters,
        delay_parameters
    )
end
