using CairoMakie, Makie, LaTeXStrings, Printf
using Gradus

include("common.jl")

function calculate_transfer_functions(m, d, angles, r0; offset = 1e-6)
    Gradus._threaded_map(angles) do angle
        x = @SVector [0.0, 1000.0, deg2rad(angle), 0.0]
        ctf = cunningham_transfer_function(m, x, d, r0, N = 500)
        mask = @. (ctf.g✶ > offset) & (ctf.g✶ < 1 - offset)
        ctf.g✶[mask], ctf.f[mask], ctf.t[mask]
    end
end

function plot_tf!(ax1, ax2, angle, X, Y, T; toff = 0.4, color = :black)
    tx = Printf.@sprintf("%0.0f", angle)

    gmax, i1max = findmax(X)
    lines!(ax1, X, Y, color = color)
    text!(
        ax1,
        gmax + 0.05,
        Y[i1max] - 0.01,
        text = L"%$(tx)^\circ",
        align = (:center, :left),
        color = color,
    )

    tt, i2max = findmax(T)
    xx = X[i2max]
    K = angle <= 35 ? 4 : 6

    lines!(ax2, X[1:i2max-K], T[1:i2max-K], color = color)
    lines!(ax2, X[i2max+K:end], T[i2max+K:end], color = color)
    text!(
        ax2,
        xx + 0.01,
        tt - toff,
        text = L"%$(tx)^\circ",
        align = (:center, :middle),
        color = color,
    )
end

m = KerrMetric(M = 1.0, a = 0.998)
d = GeometricThinDisc(0.0, 100.0, π / 2)
angles = [3, 35, 50, 65, 74, 85]

data1 = @time calculate_transfer_functions(KerrMetric(M = 1.0, a = 0.998), d, angles, 4.0)
data2 = @time calculate_transfer_functions(KerrMetric(1.0, 0.0), d, angles, 11.0)

begin
    fig = Figure(resolution = (1000, 480))
    ga = fig[1, 1] = GridLayout()
    ax1 = Axis(ga[1, 1], xlabel = L"g^\ast", ylabel = L"f", yticks = LinearTicks(4))
    ax1b = Axis(ga[1, 2], xlabel = L"g^\ast", yticks = LinearTicks(4))
    ax2 = Axis(ga[2, 1], xlabel = L"g^\ast", ylabel = L"t", yticks = LinearTicks(4))
    ax2b = Axis(ga[2, 2], xlabel = L"g^\ast", yticks = LinearTicks(4))

    palette = Iterators.Stateful(Iterators.Cycle(Makie.wong_colors()))
    for (angle, (X, Y, T)) in zip(angles, data1)
        color = popfirst!(palette)
        plot_tf!(ax1, ax2, angle, X, Y, T; color = color)
    end

    palette = Iterators.Stateful(Iterators.Cycle(Makie.wong_colors()))
    for (angle, (X, Y, T)) in zip(angles, data2)
        color = popfirst!(palette)
        plot_tf!(ax1b, ax2b, angle, X, Y, T; color = color, toff = 1.0)
    end

    Label(
        ga[1, 1, Right()],
        text = "a",
        padding = (8, 0, 180, 0),
        fontsize = 18,
        font = :bold,
    )
    Label(
        ga[1, 2, Right()],
        text = "b",
        padding = (8, 0, 180, 0),
        fontsize = 18,
        font = :bold,
    )
    Label(
        ga[2, 1, Right()],
        text = "c",
        padding = (8, 0, 180, 0),
        fontsize = 18,
        font = :bold,
    )
    Label(
        ga[2, 2, Right()],
        text = "d",
        padding = (8, 0, 180, 0),
        fontsize = 18,
        font = :bold,
    )

    hidexdecorations!(ax1, grid = false)
    hidexdecorations!(ax1b, grid = false)
    linkxaxes!(ax2, ax1)
    linkxaxes!(ax2b, ax1b)
    rowgap!(ga, 10)
    colgap!(ga, 5)

    ylims!(ax2, nothing, 1022.5)
    ylims!(ax2b, nothing, 1029.5)

    ylims!(ax1b, nothing, 0.321)

    fig
    @savefigure(fig)
end
