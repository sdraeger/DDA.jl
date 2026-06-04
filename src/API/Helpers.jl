# =============================================================================
# HELPERS
# =============================================================================

"""Write a channels × samples matrix to a temporary ASCII file for the DDA binary."""
function _write_temp_ascii(data::AbstractMatrix{<:Real})::String
    path = tempname() * ".txt"
    open(path, "w") do io
        n_ch, n_samp = size(data)
        for t in 1:n_samp
            for ch in 1:n_ch
                ch > 1 && print(io, '\t')
                @printf(io, "%.15f", Float64(data[ch, t]))
            end
            println(io)
        end
    end
    return path
end

function _resolve_labels(
    channels::AbstractVector{<:Integer},
    channel_labels::Union{Vector{String}, Nothing},
)::Vector{String}
    labels = channel_labels === nothing ? ["ch$(ch)" for ch in channels] : copy(channel_labels)
    length(labels) == length(channels) || error("Expected $(length(channels)) channel labels, got $(length(labels))")
    return labels
end

function _resolve_labels(
    file_path::AbstractString,
    channels::AbstractVector{<:Integer},
    channel_labels::Union{Vector{String}, Nothing},
)::Vector{String}
    labels = if channel_labels === nothing
        Runner._resolve_requested_channel_labels(file_path, channels; fallback_prefix="ch")
    else
        copy(channel_labels)
    end
    length(labels) == length(channels) || error("Expected $(length(channels)) channel labels, got $(length(labels))")
    return labels
end

function _pair_labels(labels::Vector{String})::Vector{String}
    pairs = String[]
    for i in 1:length(labels), j in (i + 1):length(labels)
        push!(pairs, "$(labels[i])-$(labels[j])")
    end
    return pairs
end

function _pair_channel_sets(channels::Vector{Int})::Vector{Vector{Int}}
    pairs = Vector{Vector{Int}}()
    for i in 1:length(channels), j in (i + 1):length(channels)
        push!(pairs, [channels[i], channels[j]])
    end
    return pairs
end
