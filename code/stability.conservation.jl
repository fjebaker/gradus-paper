using CairoMakie, Makie, LaTeXStrings
using Gradus

include("common.jl")

function _extract_path(sol, n_points; projection = :none, t_span = 100.0)
    mid_i = max(1, length(sol.u) ÷ 2)

    end_t = min(sol.t[mid_i] + t_span, sol.t[end])

    t_range = Gradus.Grids._inverse_grid(end_t - 18, end_t, n_points)

    r = [sol(t)[2] for t in t_range]
    θ = [sol(t)[3] for t in t_range]
    ϕ = [sol(t)[4] for t in t_range]

    if projection == :polar
        r, θ, ϕ
    else
        x = @. r * cos(ϕ) * sin(θ)
        y = @. r * sin(ϕ) * sin(θ)
        z = @. r * cos(θ)
        x, y, z
    end
end

function plot_paths_xy!(ax, sol; N = 500_000, kwargs...)
    x, y, z = _extract_path(sol, N, t_span = 16)
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
import BenchmarkTools: median
if any(i -> !isdefined(Main, i), (:t1, :t2, :t3, :t4))
    @info "Running benchmark"
    t1 = @benchmark tracegeodesics(m, x, v, 20_000.0)
    t2 = @benchmark tracegeodesics(m, x, v, 20_000.0, solver = Gradus.Feagin10())
    t3 = @benchmark tracegeodesics(m, x, v, 20_000.0, solver = Gradus.Vern6())
    t4 = @benchmark tracegeodesics(m, x, v, 20_000.0, solver = Gradus.RK4())
end

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
    translate!(axmini.scene, 0, 0, 10)
    # this needs separate translation as well, since it's drawn in the parent scene
    translate!(axmini.elements[:background], 0, 0, 9)

    hidexdecorations!(axmini)

    barplot!(axmini, [1], [median(t1.times)])
    barplot!(axmini, [2], [median(t2.times)])
    barplot!(axmini, [3], [median(t3.times)])
    barplot!(axmini, [4], [median(t4.times)])

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
    l4 = lines!(ax2, error_measure(m, x, v, solver = Gradus.RK4())...)

    leg = Legend(
        ga[1, 1],
        [l1, l2, l3, l4],
        ["Tsit5", "Feagin10", "Vern6", "RK4"],
        orientation = :horizontal,
        # labelsize = 10,
        padding = (0, 0, 0, 0),
        framevisible = false,
    )

    rowgap!(ga, 10)
    rowsize!(ga, 1, Auto(2))

    text!(ax1, (-1.3, -4.8), text = "a", font = :bold, fontsize = 20)
    text!(ax1, (11.5, 0.7), text = "b", font = :bold, fontsize = 20)
    text!(ax2, (6000, 4e-5), text = "c", font = :bold, fontsize = 20)

    resize_to_layout!(fig)
    @savefigure(fig)
    fig
end
