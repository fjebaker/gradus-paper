using CairoMakie, Makie, LaTeXStrings
using Gradus

include("common.jl")

function E_Lz(m::AbstractMetric, r)
    v = CircularOrbits.fourvelocity(m, r)
    # constrained to equitorial plane
    x = SVector(0.0, r, π / 2, 0.0)
    Gradus.E(m, x, v), Gradus.Lz(m, x, v)
end

function plot_curve!(ax, m, rs; kwargs...)
    values = E_Lz.(m, rs)
    Es = first.(values)
    Lzs = last.(values)
    lines!(ax, Lzs, Es; kwargs...)
end

m1 = KerrMetric(1.0, 0.0)
rs1 = collect(range(3.3, 50.0, 300))

m2 = KerrMetric(1.0, 0.99)
rs2 = collect(range(1.17, 50.0, 300))

begin
    fig = Figure(resolution = (500, 400))
    ax = Axis(fig[1, 1], xlabel = L"L_z", ylabel = L"E", yticks = LinearTicks(3))
    ylims!(ax, 0.71, 1.23)
    xlims!(ax, nothing, 6.5)

    plot_curve!(ax, m1, rs1, label = L"a=0.00", color = :black)
    plot_curve!(
        ax,
        m2,
        rs2,
        label = L"a=0.99",
        linestyle = :dot,
        color = :black,
        linewidth = 1.3,
    )

    axislegend(ax, position = :rb)

    fig
    @savefigure(fig)
end
