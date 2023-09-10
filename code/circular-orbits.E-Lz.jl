using CairoMakie, Makie, LaTeXStrings
using Printf
using Gradus

include("common.jl")

function E_Lz(m::AbstractMetric, r; q = 0.0)
    CircularOrbits.energy_angmom(m, SVector(r, π / 2))
end

function E_Lz(m::KerrNewmanMetric, r; q = 0.0)
    CircularOrbits.energy_angmom(m, SVector(r, π / 2); q = q)
end

function plot_curve!(ax, m, rs; q = 0.0, kwargs...)
    values = E_Lz.(m, rs; q)
    Es = first.(values)
    Lzs = last.(values)
    lines!(ax, Lzs, Es; kwargs...)
end

m1 = KerrMetric(1.0, 0.0)
rs1 = collect(range(3.3, 50.0, 300))

m2 = KerrMetric(1.0, 0.4)
rs2 = collect(range(2.5, 50.0, 300))

m3 = KerrNewmanMetric(1.0, m2.a, 1e-8)
rs3 = collect(range(2.5, 50.0, 300))
# add this one by hand
push!(rs3, 450.0)

m4 = KerrNewmanMetric(1.0, m1.a, 1e-8)
rs4 = collect(range(3.0, 50.0, 300))
# add this one by hand
push!(rs4, 450.0)

begin
    palette = Iterators.Stateful(Iterators.Cycle(Makie.wong_colors()))

    fig = Figure(resolution = (500, 400))
    layout = fig[1, 1] = GridLayout()
    ax = Axis(layout[2, 1], xlabel = L"L_z", ylabel = L"E", yticks = LinearTicks(3))
    ylims!(ax, 0.85, 1.23)
    xlims!(ax, 1, 6.5)

    color_1 = popfirst!(palette)
    _ = popfirst!(palette)
    color_2 = popfirst!(palette)
    color_3 = popfirst!(palette)

    s1 = plot_curve!(ax, m1, rs1, color = color_1)
    k1 = plot_curve!(ax, m2, rs2, linestyle = :dot, linewidth = 1.5, color = color_1)

    kp1 = plot_curve!(
        ax,
        m3,
        rs3,
        q = 0.9e8,
        color = color_2,
        linestyle = :dot,
        linewidth = 1.5,
    )
    sp1 = plot_curve!(ax, m4, rs4, q = 0.9e8, color = color_2)

    kn1 = plot_curve!(
        ax,
        m3,
        rs3,
        q = -0.9e8,
        color = color_3,
        linestyle = :dot,
        linewidth = 1.5,
    )
    sn1 = plot_curve!(ax, m4, rs4, q = -0.9e8, color = color_3)

    a_str = Printf.@sprintf("%0.1f", m2.a)

    ga = layout[1, 1] = GridLayout()
    Legend(
        ga[1, 1],
        [s1, k1],
        [L"a = 0.0", L"a = %$(a_str)"],
        orientation = :horizontal,
        framevisible = false,
        height = 10,
    )
    Legend(
        ga[2, 1],
        [sp1, s1, sn1],
        [L"qQ = 1", L"qQ = 0", L"qQ = -1"],
        orientation = :horizontal,
        framevisible = false,
    )
    colgap!(ga, 0)
    rowgap!(ga, 0)

    rowgap!(layout, 5)
    colgap!(layout, 0)

    resize_to_layout!(fig)
    @savefigure(fig)
    fig
end
