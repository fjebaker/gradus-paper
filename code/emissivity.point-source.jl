using Makie, CairoMakie, LaTeXStrings, Printf
using Gradus

include("common.jl")


function plot_reduced!(ax, x, y; kwargs...)
    N = trunc(Int, length(x) / 1000)
    I = 1:N:lastindex(x)
    lines!(ax, x[I], y[I]; kwargs...)
end

function plot_exponent!(ax, x, y; kwargs...)
    dydx = diff(y) ./ diff(x)
    x = x[2:end]
    y = y[2:end]
    em_index = @. dydx * x / y

    plot_reduced!(ax, x[1:end-1], -em_index[1:end-1]; kwargs...)
end

function calculate_emissivity(m, disc, model; kwargs...)
    prof = emissivity_profile(
        m,
        disc,
        model;
        n_samples = 70_000,
        chart = Gradus.chart_for_metric(m, closest_approach = 1.01),
        kwargs...,
    )
    isco = Gradus.isco(m)
    I = @. prof.radii >= isco
    prof.radii[I], prof.ε[I], prof.t[I]
end

function _height(model)
    v = if model.h < 2.0
        Printf.@sprintf("%.1f", model.h)
    else
        Printf.@sprintf("%.0f", model.h)
    end
    L"h=%$v r_\text{g}"
end


m = KerrMetric(1.0, 0.998)
d = ThinDisc(0.0, 1000.0)

model1 = LampPostModel(h = 1.8)
X1, Y1, T1 = @time calculate_emissivity(m, d, model1)

model2 = LampPostModel(h = 3.0)
X2, Y2, T2 = @time calculate_emissivity(m, d, model2)

model3 = LampPostModel(h = 10.0)
X3, Y3, T3 = @time calculate_emissivity(m, d, model3)

model4 = LampPostModel(h = 25.0)
X4, Y4, T4 = @time calculate_emissivity(m, d, model4)

model5 = LampPostModel(h = 100.0)
X5, Y5, T5 = @time calculate_emissivity(m, d, model5)

begin
    fig = Figure(resolution = (500, 600))

    ga = fig[1, 1] = GridLayout()
    ax1 = Axis(
        ga[1, 1],
        xscale = log10,
        yscale = log10,
        xminorgridvisible = true,
        yminorgridvisible = true,
        yticks = [1e-3, 1, 1e3],
        xticks = [1, 10, 100],
        xtickformat = values -> ["$(trunc(Int, v))" for v in values],
        ytickformat = values -> [
            value == 1 ? "1" : rich("10", superscript("$(trunc(Int,log10(value)))")) for
            value in values
        ],
        yminorticks = [1e-4, 1e-2, 1e-1, 1e1, 1e2, 1e4],
        ylabel = L"\varepsilon \,\, \text{arb.}",
    )
    ax2 = Axis(
        ga[4, 1],
        xscale = log10,
        yscale = log10,
        xminorgridvisible = true,
        xtickformat = values -> ["$(trunc(Int, v))" for v in values],
        ytickformat = values -> ["$(trunc(Int, v))" for v in values],
        xticks = [1, 10, 100],
        yminorgridvisible = true,
        yticks = [10, 100, 1000],
        yminorticks = vcat(10 .* collect(2:9), 100 .* collect(2:3)),
        ylabel = L"t",
        xlabel = L"r \, [r_\text{g}]",
    )

    ax_exponent = Axis(
        ga[3, 1],
        xscale = log10,
        xminorgridvisible = true,
        xtickformat = values -> ["$(trunc(Int, v))" for v in values],
        xticks = [1, 10, 100],
        yminorticks = [0, 2, 4, 6],
        yticks = [1, 3, 5, 7],
        yminorgridvisible = true,
        ytickformat = values -> ["    $(trunc(Int, v))" for v in values],
        ylabel = L"\alpha",
    )
    hidexdecorations!(ax_exponent, grid = false, minorgrid = false)

    linkxaxes!(ax_exponent, ax1, ax2)

    palette = _default_palette()

    c1 = popfirst!(palette)
    c2 = popfirst!(palette)
    c3 = popfirst!(palette)
    c4 = popfirst!(palette)
    _= popfirst!(palette)
    c5 = popfirst!(palette)

    l1 = plot_reduced!(ax1, X1, Y1, color = c1)
    l2 = plot_reduced!(ax1, X2, Y2, color = c2)
    l3 = plot_reduced!(ax1, X3, Y3, color = c3)
    l4 = plot_reduced!(ax1, X4, Y4, color = c4)
    l5 = plot_reduced!(ax1, X5, Y5, color = c5)

    # hlines!(ax_exponent, [3.0], color = :black, linewidth = 1.0, linestyle = :dash)
    plot_exponent!(ax_exponent, X1, Y1, color = c1)
    plot_exponent!(ax_exponent, X2, Y2, color = c2)
    plot_exponent!(ax_exponent, X3, Y3, color = c3)
    plot_exponent!(ax_exponent, X4, Y4, color = c4)
    plot_exponent!(ax_exponent, X5, Y5, color = c5)

    Legend(
        ga[2, 1],
        [l1, l2, l3, l4, l5],
        map(_height, [model1, model2, model3, model4, model5]),
        orientation = :horizontal,
        framevisible = false,
        padding = (0, 0, 0, 0),
    )

    plot_reduced!(ax2, X1, T1, color = c1)
    plot_reduced!(ax2, X2, T2, color = c2)
    plot_reduced!(ax2, X3, T3, color = c3)
    plot_reduced!(ax2, X4, T4, color = c4)
    plot_reduced!(ax2, X5, T5, color = c5)

    hidexdecorations!(ax1, grid = false, minorgrid = false)
    rowgap!(ga, 10)

    rowsize!(ga, 1, Auto(1.3))

    xlims!(ax_exponent, 0.9, 200)
    ylims!(ax2, 6, 400)
    ylims!(ax1, 1e-5, 1e4)

    ylims!(ax_exponent, -0.2, 8)

    isco = Gradus.isco(m)
    lines!(ax1, [isco, isco], [1e-5, 1e4], color = :black, linewidth = 2.0)
    lines!(ax_exponent, [isco, isco], [-1, 1e4], color = :black, linewidth = 2.0)
    lines!(ax2, [isco, isco], [6, 1e4], color = :black, linewidth = 2.0)

    Label(
        ga[2, 1, Right()],
        text = "a",
        padding = (8, 0, 155, 0),
        fontsize = 18,
        font = :bold,
    )
    Label(
        ga[3, 1, Right()],
        text = "b",
        padding = (8, 0, 100, 0),
        fontsize = 18,
        font = :bold,
    )
    Label(
        ga[4, 1, Right()],
        text = "c",
        padding = (8, 0, 100, 0),
        fontsize = 18,
        font = :bold,
    )

    @savefigure(fig)
    fig
end
