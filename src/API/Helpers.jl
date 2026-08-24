# =============================================================================
# HELPERS
# =============================================================================

"""Write a channels × samples matrix to a temporary ASCII file for the DDA binary.

The DDA ASCII reader rejects any non-numeric characters, so values must be
plain fixed-decimal floats (no `e`-notation). Formatting uses a fast
integer-scaled path for the common magnitude range and falls back to
`%.15f` elsewhere; everything is written through one buffered pass.
"""
function _write_temp_ascii(data::AbstractMatrix{<:Real})::String
    path = tempname() * ".txt"
    open(path, "w") do io
        buf = IOBuffer()
        digits = Vector{UInt8}(undef, 24)   # reusable scratch, no per-element allocs
        n_ch, n_samp = size(data)
        for t in 1:n_samp
            for ch in 1:n_ch
                ch > 1 && write(buf, UInt8('\t'))
                _print_fixed(buf, digits, Float64(data[ch, t]))
            end
            write(buf, UInt8('\n'))
        end
        write(io, take!(buf))
    end
    return path
end

# Fixed-decimal with 15 fractional digits, matching `@sprintf("%.15f", x)`
# within its documented rounding, via a single Int64 scaling. The fast path
# covers |x| < 1e3 so the scaled value stays far inside Int64 range.
const _FIXED_DECIMALS = 15
const _FIXED_SCALE_INT = 10^_FIXED_DECIMALS
const _FIXED_LIMIT = 1e3

@inline function _print_fixed(out::IOBuffer, scratch::Vector{UInt8}, x::Float64)
    if isfinite(x) && abs(x) < _FIXED_LIMIT
        r = round(Int, abs(x) * _FIXED_SCALE_INT)
        neg = signbit(x) && r != 0
        neg && write(out, UInt8('-'))
        ipart, fpart = divrem(r, _FIXED_SCALE_INT)
        # integer part (scratch holds its digits reversed)
        n = 0
        while ipart > 0
            ipart, d = divrem(ipart, 10)
            n += 1
            scratch[n] = UInt8('0' + d)
        end
        n == 0 && write(out, UInt8('0'))
        while n > 0
            write(out, scratch[n])
            n -= 1
        end
        write(out, UInt8('.'))
        # fractional part, exactly _FIXED_DECIMALS digits
        weight = _FIXED_SCALE_INT ÷ 10
        while weight > 0
            d, fpart = divrem(fpart, weight)
            write(out, UInt8('0' + d))
            weight ÷= 10
        end
    elseif isfinite(x)
        print(out, @sprintf("%.15f", x))
    elseif isnan(x)
        print(out, "nan")
    else
        print(out, x > 0 ? "inf" : "-inf")
    end
    return nothing
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

function _pair_channel_sets(channels::Vector{Int})::Vector{Vector{Int}}
    pairs = Vector{Vector{Int}}()
    for i in 1:length(channels), j in (i + 1):length(channels)
        push!(pairs, [channels[i], channels[j]])
    end
    return pairs
end
