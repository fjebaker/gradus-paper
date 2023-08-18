using CairoMakie, Makie, LaTeXStrings, Printf
using Gradus

include("common.jl")

function calculate_transfer_functions(m, d, angles, r0; offset = 1e-6)
    Gradus._threaded_map(angles) do angle
        x = @SVector [0.0, 1000.0, deg2rad(angle), 0.0]
        ctf = cunningham_transfer_function(
            m, x, d, r0, N = 500,
        )
        mask = @. (ctf.g✶ > offset) & (ctf.g✶ < 1 - offset)
        ctf.g✶[mask], ctf.f[mask]
    end
end

m = KerrMetric(M=1.0, a=0.998)
d = GeometricThinDisc(0.0, 100.0, π/2)
angles = [3, 35, 50, 65, 74, 85]

data1 = @time calculate_transfer_functions(KerrMetric(M=1.0, a=0.998), d, angles, 4.0)
data2 = @time calculate_transfer_functions(KerrMetric(1.0, 0.0), d, angles, 11.0)

begin
    fig = Figure(resolution = (500, 550))
    ga = fig[1,1] = GridLayout()
    ax = Axis(ga[1,1], xlabel = L"g^\ast", ylabel = L"f", yticks = LinearTicks(4))
    ax2 = Axis(ga[2,1], xlabel = L"g^\ast", ylabel = L"f", yticks = LinearTicks(4))

    palette = Iterators.Stateful(Iterators.Cycle(Makie.wong_colors()))
    for (angle, (X, Y)) in zip(angles, data1)
        color = popfirst!(palette)
        gmax, imax = findmax(X)
        lines!(ax, X, Y, color = color)
        tx = Printf.@sprintf("%0.0f", angle)
        text!(ax, gmax + 0.05, Y[imax] - 0.01, text = L"%$(tx)^\circ", align = (:center, :left), color = color)
    end

    palette = Iterators.Stateful(Iterators.Cycle(Makie.wong_colors()))
    for (angle, (X, Y)) in zip(angles, data2)
        color = popfirst!(palette)
        gmax, imax = findmax(X)
        lines!(ax2, X, Y, color = color)
        tx = Printf.@sprintf("%0.0f", angle)
        text!(ax2, gmax + 0.05, Y[imax] - 0.01, text = L"%$(tx)^\circ", align = (:center, :left), color = color)
    end

    Label(ga[1,1,Right()], text="a", padding = (8, 0, 180, 0), fontsize=18, font=:bold)
    Label(ga[2,1,Right()], text="b", padding = (8, 0, 180, 0), fontsize=18, font=:bold)

    hidexdecorations!(ax, grid = false)
    linkxaxes!(ax2, ax)
    rowgap!(ga, 10)

    fig
    @savefigure(fig)
end
