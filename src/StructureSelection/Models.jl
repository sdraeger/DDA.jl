"""
    make_MOD(N_MOD, DDAorder; nr_delays=2) -> Matrix{Int}

Generate the DDA structure-selection model library `MOD` used by the original
DDA scripts. `N_MOD` is the number, or list of numbers, of active monomial
terms per candidate model. `DDAorder` is the polynomial order used to construct
the `P_DDA` monomial table.
"""
function make_MOD(
    N_MOD,
    DDAorder::Integer;
    nr_delays::Integer=2,
)::Matrix{Int}
    return _structure_selection_model_space(N_MOD, DDAorder; nr_delays=nr_delays)
end
"""
    print_structure_selection([io], MOD, DDAorder; nr_delays=2, x="x")

Print the `P_DDA` monomial encoding table, `MOD` as a model-by-term checkmark
table, and each row of `MOD` as a Unicode model equation suitable for terminal
inspection.
"""
function print_structure_selection(
    io::IO,
    MOD::AbstractMatrix{<:Integer},
    DDAorder::Integer;
    nr_delays::Integer=2,
    x::AbstractString="x",
)::Nothing
    P_DDA = _p_dda(DDAorder; nr_delays=nr_delays)
    size(MOD, 2) == size(P_DDA, 1) || error(
        "MOD has $(size(MOD, 2)) columns, but P_DDA has $(size(P_DDA, 1)) monomials",
    )
    _validate_binary_MOD(MOD)

    println(io, "P_DDA")
    for idx in 1:size(P_DDA, 1)
        println(io, idx, ": ", collect(P_DDA[idx, :]))
    end

    println(io)
    _print_MOD_table(io, MOD)

    println(io)
    println(io, "Models")
    for mm in 1:size(MOD, 1)
        write_model_terminal(io, MOD, P_DDA, mm; x=x)
    end
    return nothing
end

function print_structure_selection(
    MOD::AbstractMatrix{<:Integer},
    DDAorder::Integer;
    kwargs...,
)::Nothing
    return print_structure_selection(stdout, MOD, DDAorder; kwargs...)
end

"""
    write_model_terminal([io], MOD, P_DDA, mm; x="x")

Print row `mm` of `MOD` as a Unicode equation for direct terminal inspection.
"""
function write_model_terminal(
    io::IO,
    MOD::AbstractMatrix{<:Integer},
    P_DDA::AbstractMatrix{<:Integer},
    mm::Integer;
    x::AbstractString="x",
)::Nothing
    1 <= mm <= size(MOD, 1) || error("Model row $mm out of range")
    model_terms = findall(==(1), vec(MOD[mm, :]))
    terms = P_DDA[model_terms, :]

    pieces = String[]
    for k in 1:size(terms, 1)
        monomial = _terminal_monomial(terms[k, :], x)
        push!(pieces, "a$(_subscript(k))·$monomial")
    end
    println(io, lpad(string(mm), 2), ": ẋ = ", join(pieces, " + "))
    return nothing
end

function write_model_terminal(
    MOD::AbstractMatrix{<:Integer},
    P_DDA::AbstractMatrix{<:Integer},
    mm::Integer;
    io::IO=stdout,
    kwargs...,
)::Nothing
    return write_model_terminal(io, MOD, P_DDA, mm; kwargs...)
end

"""
    write_model_LaTeX([io], MOD, SSYM, P_DDA, mm, x="x")

Print row `mm` of `MOD` as a LaTeX DDA equation using the supplied monomial
mapping table `P_DDA`.
"""
function write_model_LaTeX(
    io::IO,
    MOD::AbstractMatrix{<:Integer},
    _SSYM,
    P_DDA::AbstractMatrix{<:Integer},
    mm::Integer,
    x::AbstractString="x",
)::Nothing
    1 <= mm <= size(MOD, 1) || error("Model row $mm out of range")
    model_terms = findall(==(1), vec(MOD[mm, :]))
    terms = P_DDA[model_terms, :]
    nr_delays = maximum(P_DDA)

    @printf(io, "%2d & \\dot{%s} & = & ", mm, x)
    for k in 1:size(terms, 1)
        @printf(io, "a_%d ", k)
        for tau in 1:nr_delays
            exponent = count(==(tau), terms[k, :])
            if exponent == 1
                @printf(io, "%s_%d ", x, tau)
            elseif exponent > 1
                @printf(io, "%s_%d^%d ", x, tau, exponent)
            end
        end
        k < size(terms, 1) && print(io, "+ ")
    end
    println(io)
    return nothing
end

function write_model_LaTeX(
    MOD::AbstractMatrix{<:Integer},
    _SSYM,
    P_DDA::AbstractMatrix{<:Integer},
    mm::Integer,
    x::AbstractString="x";
    io::IO=stdout,
)::Nothing
    return write_model_LaTeX(io, MOD, _SSYM, P_DDA, mm, x)
end

function _terminal_monomial(row, x::AbstractString)::String
    parts = String[]
    nr_delays = maximum(row)
    for tau in 1:nr_delays
        exponent = count(==(tau), row)
        if exponent == 1
            push!(parts, "$x$(_subscript(tau))")
        elseif exponent > 1
            push!(parts, "$x$(_subscript(tau))$(_superscript(exponent))")
        end
    end
    return join(parts, "·")
end

function _subscript(value::Integer)::String
    return _translate_digits(value, ("₀", "₁", "₂", "₃", "₄", "₅", "₆", "₇", "₈", "₉"))
end

function _superscript(value::Integer)::String
    return _translate_digits(value, ("⁰", "¹", "²", "³", "⁴", "⁵", "⁶", "⁷", "⁸", "⁹"))
end

function _translate_digits(value::Integer, digits)::String
    value >= 0 || error("Unicode digit formatting only supports non-negative integers")
    return join(digits[Int(ch - '0') + 1] for ch in string(value))
end

function _print_MOD_table(io::IO, MOD::AbstractMatrix{<:Integer})::Nothing
    println(io, "MOD rows × P_DDA terms")
    println(io, "model | ", join(string.(1:size(MOD, 2)), " | "))
    for row_idx in 1:size(MOD, 1)
        cells = [MOD[row_idx, col_idx] == 1 ? "✓" : " " for col_idx in 1:size(MOD, 2)]
        println(io, row_idx, " | ", join(cells, " | "))
    end
    return nothing
end

function _validate_binary_MOD(MOD::AbstractMatrix{<:Integer})::Nothing
    all(value -> value == 0 || value == 1, MOD) || error("MOD must be a binary matrix")
    return nothing
end

function _structure_selection_model_space(
    N_MOD,
    DDAorder::Integer;
    nr_delays::Integer=2,
)::Matrix{Int}
    nr_delays == 2 || error("Structure-selection model generation currently supports nr_delays=2")
    DDAorder > 0 || error("DDAorder must be positive")

    P_DDA = _p_dda(DDAorder; nr_delays=nr_delays)
    model_sizes = _normalize_model_sizes(N_MOD)
    L = size(P_DDA, 1)
    mirror = _mirror_indices(P_DDA)

    rows = Vector{Vector{Int}}()
    for N in model_sizes
        1 <= N <= L || error("N_MOD entries must be between 1 and $L, got $N")
        combos = _combinations_indices(L, N)
        seen = Set{Tuple{Vararg{Int}}}()
        for combo in combos
            mirrored = Tuple(sort(mirror[combo]))
            key = min(Tuple(combo), mirrored)
            key in seen && continue
            push!(seen, key)

            row = zeros(Int, L)
            row[combo] .= 1
            push!(rows, row)
        end
    end

    isempty(rows) && error("No MOD rows generated")
    MOD = Matrix{Int}(undef, length(rows), L)
    for row_idx in eachindex(rows)
        MOD[row_idx, :] = rows[row_idx]
    end
    return MOD
end

function _p_dda(DDAorder::Integer; nr_delays::Integer=2)::Matrix{Int}
    nr_delays > 0 || error("nr_delays must be positive")
    DDAorder > 0 || error("DDAorder must be positive")
    monomials = generate_monomials(Int(nr_delays), Int(DDAorder))
    P_DDA = Matrix{Int}(undef, length(monomials), Int(DDAorder))
    for (row_idx, monomial) in enumerate(monomials)
        P_DDA[row_idx, :] = monomial
    end
    return P_DDA
end

function _normalize_model_sizes(N_MOD)::Vector{Int}
    if N_MOD isa Integer
        return [Int(N_MOD)]
    elseif N_MOD isa AbstractVector{<:Integer}
        sizes = Int.(collect(N_MOD))
        isempty(sizes) && error("N_MOD must contain at least one entry")
        return sizes
    end
    error("N_MOD must be an integer or vector/range of integers")
end

function _mirror_indices(P_DDA::AbstractMatrix{<:Integer})::Vector{Int}
    row_to_index = Dict(Tuple(P_DDA[idx, :]) => idx for idx in 1:size(P_DDA, 1))
    mirror = Vector{Int}(undef, size(P_DDA, 1))
    for idx in 1:size(P_DDA, 1)
        mirrored = _mirror_monomial(P_DDA[idx, :])
        mirror[idx] = row_to_index[Tuple(mirrored)]
    end
    return mirror
end

function _mirror_monomial(row)::Vector{Int}
    mirrored = [value == 1 ? 2 : value == 2 ? 1 : Int(value) for value in row]
    sort!(mirrored)
    return mirrored
end

function _combinations_indices(n::Integer, k::Integer)::Vector{Vector{Int}}
    result = Vector{Vector{Int}}()
    current = Int[]

    function visit(start::Int, remaining::Int)
        if remaining == 0
            push!(result, copy(current))
            return
        end
        last_start = Int(n) - remaining + 1
        for value in start:last_start
            push!(current, value)
            visit(value + 1, remaining - 1)
            pop!(current)
        end
    end

    visit(1, Int(k))
    return result
end
