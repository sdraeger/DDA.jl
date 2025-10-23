using Dates

struct TimeRange
    start::Float64
    stop::Float64
end

struct PreprocessingOptions
    detrending::Union{String,Nothing}
    highpass::Union{Float64,Nothing}
    lowpass::Union{Float64,Nothing}
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

struct ScaleParameters
    scale_min::Float64
    scale_max::Float64
    scale_num::Int
end

struct DDARequest
    file_path::String
    channels::Union{Vector{Int},Nothing}
    time_range::TimeRange
    preprocessing_options::PreprocessingOptions
    algorithm_selection::AlgorithmSelection
    window_parameters::WindowParameters
    scale_parameters::ScaleParameters
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
    scale_parameters::ScaleParameters
    created_at::DateTime
end

function DDAResult(
    id::String,
    file_path::String,
    channels::Vector{String},
    q_matrix::Matrix{Float64},
    window_parameters::WindowParameters,
    scale_parameters::ScaleParameters
)
    DDAResult(
        id,
        file_path,
        channels,
        q_matrix,
        nothing,
        nothing,
        window_parameters,
        scale_parameters,
        now()
    )
end
