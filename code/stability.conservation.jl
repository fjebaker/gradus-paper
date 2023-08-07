using CairoMakie, Makie, LaTeXStrings
using Gradus

include("common.jl")

function plot_paths_xy!(ax, sol; N = 200_000, kwargs...)
    x, y, z = Gradus._extract_path(sol, N, t_span = 16)
    lines!(ax, x, y; kwargs...)
end

function error_measure(m, x, v; kwargs...)
    sol = tracegeodesics(m, x, v, 20_000.0; kwargs...)

    points = unpack_solution_full(m, sol)
    dotprod = [Gradus.dotproduct(m, p.x, p.v, p.v) for p in points]
    rs = [p.x[2] for p in points]
    rs, abs.(dotprod)
end

m = KerrMetric(a = 0.998)
x = SVector(0.0, 10_000.0, π / 2, 0.0)
v = map_impact_parameters(m, x, 5.0, 0.0)

sol = tracegeodesics(m, x, v, 20_000.0)

using BenchmarkTools
t1 = @benchmark tracegeodesics(m, x, v, 20_000.0)
t2 = @benchmark tracegeodesics(m, x, v, 20_000.0, solver = Gradus.Feagin10())
t3 = @benchmark tracegeodesics(m, x, v, 20_000.0, solver = Gradus.Vern6())

begin
    fig = Figure(resolution = (480, 550))
    ga = fig[1, 1] = GridLayout()

    ax1 = Axis(
        ga[2, 1],
        aspect = DataAspect(),
        xticks = LinearTicks(4),
        yticks = LinearTicks(4),
        xlabel = L"x\, [r_\text{g}]",
        ylabel = L"y\, [r_\text{g}]",
    )
    ylims!(ax1, -5, 2)
    xlims!(ax1, -2, 13)

    plot_paths_xy!(ax1, sol)

    axmini = Axis(
        ga[2, 1],
        width = Relative(0.35),
        height = Relative(0.5),
        halign = 0.8,
        valign = 0.8,
        xgridvisible = false,
        ygridvisible = false,
        ytickformat = values -> ["$(trunc(Int, v / 1e3)) μs" for v in values],
    )

    hidexdecorations!(axmini)

    barplot!(axmini, [1], [mean(t1.times)])
    barplot!(axmini, [2], [mean(t2.times)])
    barplot!(axmini, [3], [mean(t3.times)])

    ax2 = Axis(
        ga[3, 1],
        yscale = log10,
        xscale = log10,
        ylabel = L"\text{abs} \left(\, g_{\mu\nu} v^\mu v^\nu  \,  \right)",
        xlabel = L"r\, [r_\text{g}]",
        xtickformat = values -> ["$(trunc(Int, v))" for v in values],
    )

    l1 = lines!(ax2, error_measure(m, x, v)...)
    l2 = lines!(ax2, error_measure(m, x, v, solver = Gradus.Feagin10())...)
    l3 = lines!(ax2, error_measure(m, x, v, solver = Gradus.Vern6())...)

    leg = Legend(
        ga[1, 1],
        [l1, l2, l3],
        ["Tsit5", "Feagin10", "Vern6"],
        orientation = :horizontal,
        labelsize = 10,
        padding = (0, 0, 0, 0),
        framevisible = false,
    )

    rowgap!(ga, 10)
    rowsize!(ga, 1, Auto(2))

    text!(ax1, (-1.3, -4.8), text = "a", font = :bold, fontsize = 20)
    text!(ax1, (11.5, 0.7), text = "b", font = :bold, fontsize = 20)
    text!(ax2, (6000, 4e-7), text = "c", font = :bold, fontsize = 20)

    # resize_to_layout!(fig)
    fig
    @savefigure(fig)
end
