function _fallback_channel_label(channel::Integer, fallback_prefix::String)::String
    return string(fallback_prefix, Int(channel))
end

function _sanitize_channel_label(label::AbstractString)::String
    cleaned = replace(String(label), '\0' => ' ')
    cleaned = strip(cleaned)
    isempty(cleaned) && return ""
    return strip(cleaned, ['"', '\''])
end

function _read_edf_channel_labels(file_path::AbstractString)::Union{Vector{String}, Nothing}
    open(file_path, "r") do io
        fixed_header = read(io, 256)
        length(fixed_header) == 256 || return nothing

        signal_count = tryparse(Int, strip(String(fixed_header[253:256])))
        signal_count === nothing && return nothing
        signal_count > 0 || return nothing

        labels = String[]
        for _ in 1:signal_count
            field = read(io, 16)
            length(field) == 16 || return nothing
            push!(labels, _sanitize_channel_label(String(field)))
        end
        return labels
    end
end

function _split_ascii_fields(line::AbstractString)::Vector{String}
    stripped = strip(replace(line, '\ufeff' => ' '))
    isempty(stripped) && return String[]

    if occursin('\t', stripped)
        parts = split(stripped, '\t'; keepempty=true)
    elseif occursin(',', stripped)
        parts = split(stripped, ','; keepempty=true)
    else
        parts = split(stripped)
    end

    return [_sanitize_channel_label(part) for part in parts]
end

function _is_numeric_field(field::AbstractString)::Bool
    stripped = strip(field)
    isempty(stripped) && return false
    return tryparse(Float64, stripped) !== nothing
end

function _read_ascii_channel_labels(file_path::AbstractString)::Union{Vector{String}, Nothing}
    return open(file_path, "r") do io
        for line in eachline(io)
            stripped = strip(replace(line, '\ufeff' => ' '))
            (isempty(stripped) || startswith(stripped, '#')) && continue

            fields = _split_ascii_fields(stripped)
            isempty(fields) && continue

            all(_is_numeric_field, fields) && return nothing
            return fields
        end

        return nothing
    end
end

function _infer_input_channel_labels(file_path::AbstractString)::Union{Vector{String}, Nothing}
    isfile(file_path) || return nothing

    try
        ext = lowercase(splitext(String(file_path))[2])
        labels = ext == ".edf" ? _read_edf_channel_labels(file_path) : _read_ascii_channel_labels(file_path)
        labels === nothing && return nothing
        any(!isempty, labels) || return nothing
        return labels
    catch
        return nothing
    end
end

function _resolve_requested_channel_labels(
    file_path::AbstractString,
    channels::AbstractVector{<:Integer};
    fallback_prefix::String="Channel ",
)::Vector{String}
    inferred = _infer_input_channel_labels(file_path)
    resolved = String[]

    for channel in channels
        idx = Int(channel)
        if inferred !== nothing && idx <= length(inferred)
            label = _sanitize_channel_label(inferred[idx])
            if !isempty(label)
                push!(resolved, label)
                continue
            end
        end
        push!(resolved, _fallback_channel_label(idx, fallback_prefix))
    end

    return resolved
end
