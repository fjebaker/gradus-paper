using CairoMakie, Makie, LaTeXStrings
using Printf
using Gradus

include("common.jl")

function E_Lz(m::AbstractMetric, r; q = 0.0)
    CircularOrbits.energy_angmom(m, SVector(r, π / 2))
end

function plot_curve!(ax, m, rs; q = 0.0, kwargs...)
    values = E_Lz.(m, rs; q)
    Es = first.(values)
    Lzs = last.(values)
    lines!(ax, Lzs, Es; kwargs...)
end

m0 = KerrMetric(1.0, 0.0)
rs0 = collect(range(3.3, 50.0, 300))

m1 = KerrMetric(1.0, 0.4)
rs1 = collect(range(2.5, 50.0, 300))

m2 = KerrMetric(1.0, 0.8)
rs2 = collect(range(1.96, 50.0, 300))

m3 = KerrMetric(1.0, 0.998)
rs3 = collect(Gradus.Grids._geometric_grid(1.08, 50.0, 300))

begin
    palette = _default_palette()

    fig = Figure(resolution = (450, 300))
    layout = fig[1, 1] = GridLayout()
    ax = Axis(layout[2, 1], xlabel = L"L_z", ylabel = L"E", yticks = LinearTicks(3))
    ylims!(ax, 0.65, 1.15)
    xlims!(ax, 1, 6.5)

    s0 = plot_curve!(ax, m0, rs0, color = popfirst!(palette))
    k1 = plot_curve!(ax, m1, rs1, linewidth = 1.5, color = popfirst!(palette))
    k2 = plot_curve!(ax, m2, rs2, linewidth = 1.5, color = popfirst!(palette))
    k3 = plot_curve!(ax, m3, rs3, linewidth = 1.5, color = popfirst!(palette))

    a1_str = Printf.@sprintf("%0.1f", m1.a)
    a2_str = Printf.@sprintf("%0.1f", m2.a)
    a3_str = Printf.@sprintf("%0.1f", m3.a)

    ga = layout[1, 1] = GridLayout()
    Legend(
        ga[1, 1],
        [k3, k2, k1, s0],
        [L"a = %$(a)" for a in reverse(("0.0", a1_str, a2_str, a3_str))],
        orientation = :horizontal,
        framevisible = false,
        height = 10,
    )
    colgap!(ga, 0)
    rowgap!(ga, 0)

    rowgap!(layout, 5)
    colgap!(layout, 0)

    resize_to_layout!(fig)
    @savefigure(fig)
    fig
end
