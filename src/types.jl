using Dates

struct TimeRange
    start::Float64
    stop::Float64
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
    delay_min::Float64
    delay_max::Float64
    delay_num::Int
end

struct DDARequest
    file_path::String
    channels::Union{Vector{Int},Nothing}
    time_range::TimeRange
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
    id::String
    file_path::String
    channels::Vector{String}
    q_matrix::Matrix{Float64}
    variant_results::Union{Vector{VariantResult},Nothing}
    raw_output::Union{String,Nothing}
    window_parameters::WindowParameters
    delay_parameters::DelayParameters
end

function DDAResult(
    id::String,
    file_path::String,
    channels::Vector{String},
    q_matrix::Matrix{Float64},
    window_parameters::WindowParameters,
    delay_parameters::DelayParameters
)
    DDAResult(
        id,
        file_path,
        channels,
        q_matrix,
        nothing,
        nothing,
        window_parameters,
        delay_parameters
    )
end
