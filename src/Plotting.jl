"""Plotting functions for DDA results (lazy Plots.jl loading)."""
module Plotting

using ..Results: STResult, CTResult, DEResult, n_windows
using ..ModelEncoding: generate_monomials, monomial_to_text
using ..StructureSelection
using ..StructureSelection: StructureSelectionRun
using ..OptionalDeps

# =============================================================================
# Lazy Plots.jl loading
# =============================================================================

_ensure_plots()::Module = OptionalDeps.require(:Plots)

_entity_labels(result::STResult) = result.channel_labels
_entity_labels(result::CTResult) = result.pair_labels

function _plot_x_axis(result; use_time::Bool, sfreq::Union{Float64,Nothing})
    if use_time && !isempty(result.t)
        return result.t, "t"
    elseif use_time && sfreq !== nothing
        return result.window_starts ./ sfreq, "Time (s)"
    end
    return 1:n_windows(result), "Window"
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
    P = _ensure_plots()
    coeffs = result.coefficients
    nc = size(coeffs, 1)
    nk = size(coeffs, 3)

    ch_idx = something(channels, collect(1:nc))
    k_idx = something(coeff_indices, collect(1:nk))

    labels = _entity_labels(result)
    x, xlabel = _plot_x_axis(result; use_time=use_time, sfreq=sfreq)

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
    P = _ensure_plots()
    nc = size(result.errors, 1)
    ch_idx = something(channels, collect(1:nc))
    labels = _entity_labels(result)
    x, xlabel = _plot_x_axis(result; use_time=use_time, sfreq=sfreq)

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
    P = _ensure_plots()
    data_2d = result.coefficients[:, :, coeff_index]
    labels = _entity_labels(result)
    x, xlabel = _plot_x_axis(result; use_time=use_time, sfreq=sfreq)

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
    P = _ensure_plots()
    x, xlabel = _plot_x_axis(result; use_time=use_time, sfreq=sfreq)

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
    P = _ensure_plots()
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

# =============================================================================
# plot_structure_selection
# =============================================================================

"""
    plot_structure_selection(run; mode=:all, channels=nothing, models=nothing,
                             metric=:mean_error, figsize=(800,600))

Plot cached two-delay structure-selection results. `mode=:all` shows the
winning model at each `(tau1, tau2)` pair after aggregating over all windows.
`mode=:time` shows the winning delay pair at each raw DDA window coordinate.
Model colors are the corresponding row numbers in `run.MOD`.
"""
function plot_structure_selection(
    run::StructureSelectionRun;
    mode::Symbol=:all,
    channels=nothing,
    models=nothing,
    metric::Symbol=:mean_error,
    figsize::Tuple{Int,Int}=(800, 600),
)
    _ensure_plots()
    data = StructureSelection._structure_selection_plot_data(
        run;
        mode=mode,
        channels=channels,
        models=models,
        metric=metric,
    )
    return Base.invokelatest(_plot_structure_selection_impl, data; figsize=figsize)
end

function _plot_structure_selection_impl(data; figsize)
    P = _ensure_plots()
    models = sort(unique(filter(>(0), vec(data.model_numbers))))
    model_indices = Dict(model => idx for (idx, model) in enumerate(models))
    palette_colors = [
        "#35618F",
        "#C65D3B",
        "#4F8A70",
        "#9B6AA6",
        "#C08A35",
        "#4F8C99",
        "#A5545D",
        "#6F7782",
    ]
    palette = length(models) <= length(palette_colors) ?
              Base.invokelatest(P.palette, palette_colors[1:length(models)]) :
              Base.invokelatest(P.palette, :tab20, length(models))
    colors = Base.invokelatest(P.cgrad, palette; categorical=true)
    clims = (0.5, length(models) + 0.5)
    style = (
        background_color=:white,
        foreground_color_axis="#333333",
        foreground_color_border="#333333",
        foreground_color_text="#262626",
        fontfamily="sans-serif",
        framestyle=:box,
        guidefontsize=11,
        tickfontsize=9,
        legendfontsize=9,
        legend_background_color=:transparent,
        legend_foreground_color=:transparent,
        tick_direction=:out,
        dpi=300,
    )

    if data.mode == :all
        values = fill(NaN, size(data.model_numbers))
        for idx in eachindex(values)
            model = data.model_numbers[idx]
            model > 0 && (values[idx] = model_indices[model])
        end
        plot = Base.invokelatest(
            P.heatmap,
            data.tau1,
            data.tau2,
            values;
            style...,
            xlabel="Delay τ₁",
            ylabel="Delay τ₂",
            xticks=_structure_selection_ticks(data.tau1),
            yticks=_structure_selection_ticks(data.tau2),
            color=colors,
            clims=clims,
            colorbar=false,
            legend=:outerright,
            aspect_ratio=:equal,
            grid=false,
            size=figsize,
        )
        for model in models
            Base.invokelatest(
                P.scatter!,
                plot,
                [NaN],
                [NaN];
                color=palette[model_indices[model]],
                markershape=:rect,
                markersize=6,
                markerstrokewidth=0,
                label="Model $model",
            )
        end
        return plot
    end

    delays = sort(unique(vcat(data.tau1, data.tau2)))
    padding = max(1.0, 0.05 * (maximum(delays) - minimum(delays)))
    ylims = (minimum(delays) - padding, maximum(delays) + padding)
    p1 = Base.invokelatest(
        P.plot;
        style...,
        ylabel="Delay τ₁",
        xticks=false,
        yticks=_structure_selection_ticks(delays),
        ylims=ylims,
        grid=:y,
        gridalpha=0.16,
        gridlinewidth=0.6,
        legend=false,
    )
    p2 = Base.invokelatest(
        P.plot;
        style...,
        xlabel="Window coordinate, T",
        ylabel="Delay τ₂",
        yticks=_structure_selection_ticks(delays),
        ylims=ylims,
        grid=:y,
        gridalpha=0.16,
        gridlinewidth=0.6,
        legend=:outerright,
    )
    Base.invokelatest(
        P.plot!,
        p1,
        data.T,
        data.tau1;
        color="#AEB4BA",
        linewidth=1.0,
        seriestype=:steppost,
        label="",
    )
    Base.invokelatest(
        P.plot!,
        p2,
        data.T,
        data.tau2;
        color="#AEB4BA",
        linewidth=1.0,
        seriestype=:steppost,
        label="",
    )
    for model in models
        indices = findall(==(model), data.model_numbers)
        color = palette[model_indices[model]]
        Base.invokelatest(
            P.scatter!,
            p1,
            data.T[indices],
            data.tau1[indices];
            color=color,
            markersize=4,
            markerstrokecolor="#333333",
            markerstrokewidth=0.5,
            label="",
        )
        Base.invokelatest(
            P.scatter!,
            p2,
            data.T[indices],
            data.tau2[indices];
            color=color,
            markersize=4,
            markerstrokecolor="#333333",
            markerstrokewidth=0.5,
            label="Model $model",
        )
    end
    return Base.invokelatest(P.plot, p1, p2; layout=(2, 1), link=:x, size=figsize)
end

function _structure_selection_ticks(values; maximum_ticks::Int=8)
    length(values) <= maximum_ticks && return values
    indices = unique(round.(Int, range(1, length(values), length=maximum_ticks)))
    return values[indices]
end

end # module Plotting
