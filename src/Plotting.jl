"""Plotting functions for DDA results (lazy Plots.jl loading)."""
module Plotting

using ..Results
using ..ModelEncoding

export plot_coefficients, plot_heatmap, plot_errors, plot_ergodicity, plot_model

# =============================================================================
# Lazy Plots.jl loading
# =============================================================================

const _plots_loaded = Ref(false)

function _ensure_plots()
    if !_plots_loaded[]
        try
            @eval import Plots
            _plots_loaded[] = true
        catch
            error(
                "Plots.jl is required for plotting functions. " *
                "Install with: using Pkg; Pkg.add(\"Plots\")"
            )
        end
    end
end

# =============================================================================
# plot_coefficients
# =============================================================================

"""
    plot_coefficients(result; coeff_indices=nothing, channels=nothing,
                      use_time=false, sfreq=nothing, figsize=(800,400))

Plot DDA coefficient time series.

Requires `Plots.jl` to be installed.

# Arguments
- `result::Union{STResult, CTResult}`: Analysis result.
- `coeff_indices`: Which coefficients to plot (1-based). Default: all.
- `channels`: Which channels to plot (1-based). Default: all.
- `use_time::Bool=false`: Use `result.t` when available. Falls back to `window_starts / sfreq`.
- `sfreq`: Sampling frequency for fallback time conversion.
- `figsize`: Figure size as (width, height).
"""
function plot_coefficients(
    result::Union{STResult,CTResult};
    coeff_indices::Union{Vector{Int},Nothing}=nothing,
    channels::Union{Vector{Int},Nothing}=nothing,
    use_time::Bool=false,
    sfreq::Union{Float64,Nothing}=nothing,
    figsize::Tuple{Int,Int}=(800, 400),
)
    _ensure_plots()
    Base.invokelatest(_plot_coefficients_impl, result;
        coeff_indices=coeff_indices, channels=channels,
        use_time=use_time, sfreq=sfreq, figsize=figsize)
end

function _plot_coefficients_impl(
    result::Union{STResult,CTResult};
    coeff_indices, channels, use_time, sfreq, figsize,
)
    P = @eval Plots
    coeffs = result.coefficients
    nc = size(coeffs, 1)
    nk = size(coeffs, 3)

    ch_idx = something(channels, collect(1:nc))
    k_idx = something(coeff_indices, collect(1:nk))

    labels = result isa STResult ? result.channel_labels : result.pair_labels

    if use_time && !isempty(result.t)
        x = result.t
        xlabel = "t"
    elseif use_time && sfreq !== nothing
        x = result.window_starts ./ sfreq
        xlabel = "Time (s)"
    else
        x = 1:n_windows(result)
        xlabel = "Window"
    end

    plots = []
    for ki in k_idx
        p = Base.invokelatest(P.plot;
            title="a_$ki", xlabel=xlabel, ylabel="a_$ki",
            size=figsize, legend=:outerright)
        for ci in ch_idx
            Base.invokelatest(P.plot!, p, x, coeffs[ci, :, ki];
                label=labels[ci])
        end
        push!(plots, p)
    end

    if length(plots) == 1
        return plots[1]
    else
        return Base.invokelatest(P.plot, plots...; layout=(length(plots), 1),
            size=(figsize[1], figsize[2] * length(plots)))
    end
end

# =============================================================================
# plot_errors
# =============================================================================

"""
    plot_errors(result; channels=nothing, use_time=false, sfreq=nothing, figsize=(800,400))

Plot reconstruction error time series.
"""
function plot_errors(
    result::Union{STResult,CTResult};
    channels::Union{Vector{Int},Nothing}=nothing,
    use_time::Bool=false,
    sfreq::Union{Float64,Nothing}=nothing,
    figsize::Tuple{Int,Int}=(800, 400),
)
    _ensure_plots()
    Base.invokelatest(_plot_errors_impl, result;
        channels=channels, use_time=use_time, sfreq=sfreq, figsize=figsize)
end

function _plot_errors_impl(
    result::Union{STResult,CTResult};
    channels, use_time, sfreq, figsize,
)
    P = @eval Plots
    nc = size(result.errors, 1)
    ch_idx = something(channels, collect(1:nc))
    labels = result isa STResult ? result.channel_labels : result.pair_labels

    if use_time && !isempty(result.t)
        x = result.t
        xlabel = "t"
    elseif use_time && sfreq !== nothing
        x = result.window_starts ./ sfreq
        xlabel = "Time (s)"
    else
        x = 1:n_windows(result)
        xlabel = "Window"
    end

    p = Base.invokelatest(P.plot;
        title="Reconstruction Error", xlabel=xlabel, ylabel="Error",
        size=figsize, legend=:outerright)
    for ci in ch_idx
        Base.invokelatest(P.plot!, p, x, result.errors[ci, :];
            label=labels[ci])
    end
    return p
end

# =============================================================================
# plot_heatmap
# =============================================================================

"""
    plot_heatmap(result; coeff_index=1, use_time=false, sfreq=nothing,
                 cmap=:RdBu, figsize=(800,400))

Plot a heatmap of one coefficient across channels and windows.
"""
function plot_heatmap(
    result::Union{STResult,CTResult};
    coeff_index::Int=1,
    use_time::Bool=false,
    sfreq::Union{Float64,Nothing}=nothing,
    cmap::Symbol=:RdBu,
    vmin::Union{Float64,Nothing}=nothing,
    vmax::Union{Float64,Nothing}=nothing,
    figsize::Tuple{Int,Int}=(800, 400),
)
    _ensure_plots()
    Base.invokelatest(_plot_heatmap_impl, result;
        coeff_index=coeff_index, use_time=use_time, sfreq=sfreq,
        cmap=cmap, vmin=vmin, vmax=vmax, figsize=figsize)
end

function _plot_heatmap_impl(
    result::Union{STResult,CTResult};
    coeff_index, use_time, sfreq, cmap, vmin, vmax, figsize,
)
    P = @eval Plots
    data_2d = result.coefficients[:, :, coeff_index]
    labels = result isa STResult ? result.channel_labels : result.pair_labels

    if use_time && !isempty(result.t)
        x = result.t
        xlabel = "t"
    elseif use_time && sfreq !== nothing
        x = result.window_starts ./ sfreq
        xlabel = "Time (s)"
    else
        x = 1:n_windows(result)
        xlabel = "Window"
    end

    abs_max = maximum(abs, data_2d)
    clims = (something(vmin, -abs_max), something(vmax, abs_max))

    p = Base.invokelatest(P.heatmap, x, 1:size(data_2d, 1), data_2d;
        title="a_$coeff_index",
        xlabel=xlabel, ylabel="Channel",
        yticks=(1:length(labels), labels),
        color=cmap, clims=clims,
        size=figsize)
    return p
end

# =============================================================================
# plot_ergodicity
# =============================================================================

"""
    plot_ergodicity(result; use_time=false, sfreq=nothing, figsize=(800,400))

Plot ergodicity measure from a DEResult.
"""
function plot_ergodicity(
    result::DEResult;
    use_time::Bool=false,
    sfreq::Union{Float64,Nothing}=nothing,
    figsize::Tuple{Int,Int}=(800, 400),
)
    _ensure_plots()
    Base.invokelatest(_plot_ergodicity_impl, result;
        use_time=use_time, sfreq=sfreq, figsize=figsize)
end

function _plot_ergodicity_impl(
    result::DEResult;
    use_time, sfreq, figsize,
)
    P = @eval Plots

    if use_time && !isempty(result.t)
        x = result.t
        xlabel = "t"
    elseif use_time && sfreq !== nothing
        x = result.window_starts ./ sfreq
        xlabel = "Time (s)"
    else
        x = 1:n_windows(result)
        xlabel = "Window"
    end

    p = Base.invokelatest(P.plot, x, result.ergodicity;
        title="Dynamical Ergodicity",
        xlabel=xlabel, ylabel="Ergodicity",
        linewidth=1.5, color=:steelblue,
        legend=false, size=figsize)
    return p
end

# =============================================================================
# plot_model
# =============================================================================

"""
    plot_model(model_encoding; num_delays=2, polynomial_order=4, figsize=(600,400))

Visualize the DDA model space grid with selected terms highlighted.
"""
function plot_model(
    model_encoding::Vector{Int};
    num_delays::Int=2,
    polynomial_order::Int=4,
    figsize::Tuple{Int,Int}=(600, 400),
)
    _ensure_plots()
    Base.invokelatest(_plot_model_impl, model_encoding;
        num_delays=num_delays, polynomial_order=polynomial_order, figsize=figsize)
end

function _plot_model_impl(
    model_encoding::Vector{Int};
    num_delays, polynomial_order, figsize,
)
    P = @eval Plots
    monomials = generate_monomials(num_delays, polynomial_order)
    selected = Set(model_encoding)

    # Group monomials by effective degree
    degrees = Dict{Int, Vector{Tuple{Int,String}}}()
    for (idx, mono) in enumerate(monomials)
        deg = count(x -> x != 0, mono)
        term = monomial_to_text(mono)
        if !haskey(degrees, deg)
            degrees[deg] = Tuple{Int,String}[]
        end
        push!(degrees[deg], (idx, term))
    end

    max_deg = maximum(keys(degrees))
    max_cols = maximum(length(v) for v in values(degrees))

    p = Base.invokelatest(P.plot;
        title="Model Space (delays=$num_delays, order=$polynomial_order)",
        xlim=(0, max_cols + 1), ylim=(0, max_deg + 1),
        size=figsize, legend=false,
        xticks=nothing, grid=false, framestyle=:none)

    for deg in 1:max_deg
        haskey(degrees, deg) || continue
        for (col, (idx, term)) in enumerate(degrees[deg])
            color = idx in selected ? :green : :lightgray
            Base.invokelatest(P.plot!, p,
                [col - 0.4, col + 0.4, col + 0.4, col - 0.4, col - 0.4],
                [deg - 0.4, deg - 0.4, deg + 0.4, deg + 0.4, deg - 0.4];
                fill=true, fillcolor=color, linecolor=:gray, fillalpha=0.7)
            Base.invokelatest(P.annotate!, p, col, deg,
                Base.invokelatest(P.text, "[$idx] $term", 7))
        end
    end

    Base.invokelatest(P.ylabel!, p, "Degree")
    return p
end

end # module Plotting
