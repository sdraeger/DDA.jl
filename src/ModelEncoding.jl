"""Model space encoding utilities for DDA.

Provides functions to enumerate monomials, decode model encodings to
human-readable equations, and visualize the model space.
"""
module ModelEncoding

using Printf

export generate_monomials, monomial_to_text, monomial_to_latex
export decode_model_encoding, visualize_model_space

"""
    generate_monomials(num_delays, polynomial_order) -> Vector{NTuple{N,Int}}

Generate all non-decreasing tuples of length `polynomial_order` with values
in `{0, 1, ..., num_delays}`, excluding the all-zeros tuple.

Each tuple represents a monomial where value 0 = "no factor" (padding) and
value k > 0 = x(t - τ_k).

# Examples
```julia
generate_monomials(2, 2)
# [(0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]
```
"""
function generate_monomials(num_delays::Int, polynomial_order::Int)::Vector{Vector{Int}}
    results = Vector{Vector{Int}}()
    _non_decreasing_tuples!(results, Int[], num_delays, polynomial_order, 0)
    # Remove the all-zeros tuple
    filter!(t -> any(x -> x != 0, t), results)
    return results
end

function _non_decreasing_tuples!(
    results::Vector{Vector{Int}},
    current::Vector{Int},
    num_delays::Int,
    remaining::Int,
    min_val::Int,
)
    if remaining == 0
        push!(results, copy(current))
        return
    end
    for v in min_val:num_delays
        push!(current, v)
        _non_decreasing_tuples!(results, current, num_delays, remaining - 1, v)
        pop!(current)
    end
end

"""
    monomial_to_text(monomial; tau_values=nothing) -> String

Convert a monomial tuple to a plain-text representation.

# Examples
```julia
monomial_to_text([0, 1])       # "x_1"
monomial_to_text([1, 1])       # "x_1^2"
monomial_to_text([1, 2])       # "x_1 * x_2"
monomial_to_text([0, 1]; tau_values=[7, 10])  # "x(t-7)"
```
"""
function monomial_to_text(monomial::Vector{Int}; tau_values::Union{Vector{Int},Nothing}=nothing)::String
    factors = filter(x -> x != 0, monomial)
    isempty(factors) && return "1"

    counts = Dict{Int,Int}()
    for f in factors
        counts[f] = get(counts, f, 0) + 1
    end

    parts = String[]
    for k in sort(collect(keys(counts)))
        exp = counts[k]
        if tau_values !== nothing && k <= length(tau_values)
            base = "x(t-$(tau_values[k]))"
        else
            base = "x_$k"
        end
        if exp == 1
            push!(parts, base)
        else
            push!(parts, "$(base)^$exp")
        end
    end
    return join(parts, " * ")
end

"""
    monomial_to_latex(monomial; tau_values=nothing) -> String

Convert a monomial tuple to a LaTeX representation.
"""
function monomial_to_latex(monomial::Vector{Int}; tau_values::Union{Vector{Int},Nothing}=nothing)::String
    factors = filter(x -> x != 0, monomial)
    isempty(factors) && return "1"

    counts = Dict{Int,Int}()
    for f in factors
        counts[f] = get(counts, f, 0) + 1
    end

    parts = String[]
    for k in sort(collect(keys(counts)))
        exp = counts[k]
        if tau_values !== nothing && k <= length(tau_values)
            base = "x(t-$(tau_values[k]))"
        else
            base = "x_{$k}"
        end
        if exp == 1
            push!(parts, base)
        else
            push!(parts, "$(base)^{$exp}")
        end
    end
    return join(parts, " \\cdot ")
end

"""
    decode_model_encoding(model_encoding; num_delays=2, polynomial_order=4,
                          tau_values=nothing, format="text") -> String

Decode a model encoding (e.g. `[1, 2, 10]`) into a human-readable DDE equation.

# Arguments
- `model_encoding::Vector{Int}`: 1-based indices into the monomial space.
- `num_delays::Int=2`: Number of delay variables.
- `polynomial_order::Int=4`: Maximum polynomial order.
- `tau_values::Union{Vector{Int},Nothing}=nothing`: Actual delay values for display.
- `format::String="text"`: Output format, `"text"` or `"latex"`.

# Examples
```julia
decode_model_encoding([1, 2, 10]; num_delays=2, polynomial_order=4)
# "dx/dt = a_1 x_1 + a_2 x_2 + a_3 x_1^4"
```
"""
function decode_model_encoding(
    model_encoding::Vector{Int};
    num_delays::Int=2,
    polynomial_order::Int=4,
    tau_values::Union{Vector{Int},Nothing}=nothing,
    format::String="text",
)::String
    monomials = generate_monomials(num_delays, polynomial_order)
    n_mono = length(monomials)

    for idx in model_encoding
        (1 <= idx <= n_mono) || error("Model index $idx out of range [1, $n_mono]")
    end

    to_str = format == "latex" ? monomial_to_latex : monomial_to_text

    terms = String[]
    for (i, idx) in enumerate(model_encoding)
        mono = monomials[idx]
        term_str = to_str(mono; tau_values=tau_values)
        if format == "latex"
            push!(terms, "a_{$i} $(term_str)")
        else
            push!(terms, "a_$i $(term_str)")
        end
    end

    lhs = format == "latex" ? "\\frac{dx}{dt}" : "dx/dt"
    return "$lhs = " * join(terms, " + ")
end

"""
    visualize_model_space(num_delays, polynomial_order; tau_values=nothing,
                          highlight_encoding=nothing) -> String

Generate a text table of the full model space with optional highlighting
of selected model terms.

# Examples
```julia
println(visualize_model_space(2, 4; highlight_encoding=[1, 2, 10]))
```
"""
function visualize_model_space(
    num_delays::Int,
    polynomial_order::Int;
    tau_values::Union{Vector{Int},Nothing}=nothing,
    highlight_encoding::Union{Vector{Int},Nothing}=nothing,
)::String
    monomials = generate_monomials(num_delays, polynomial_order)
    n = length(monomials)

    lines = String[]
    push!(lines, "Model Space: num_delays=$num_delays, polynomial_order=$polynomial_order")
    push!(lines, "Total monomials: $n")
    push!(lines, "")
    push!(lines, @sprintf("%-6s %-20s %-30s", "Index", "Encoding", "Term"))
    push!(lines, repeat("-", 60))

    highlighted = highlight_encoding === nothing ? Set{Int}() : Set(highlight_encoding)

    for (i, mono) in enumerate(monomials)
        term = monomial_to_text(mono; tau_values=tau_values)
        marker = i in highlighted ? " *" : ""
        push!(lines, @sprintf("%-6d %-20s %-30s%s", i, string(mono), term, marker))
    end

    if highlight_encoding !== nothing
        push!(lines, "")
        push!(lines, "Selected terms (marked with *):")
        eq = decode_model_encoding(highlight_encoding;
            num_delays=num_delays, polynomial_order=polynomial_order,
            tau_values=tau_values, format="text")
        push!(lines, eq)
    end

    return join(lines, "\n")
end

end # module ModelEncoding
