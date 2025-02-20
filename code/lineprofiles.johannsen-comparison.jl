using Gradus
using Printf

include("common.jl")

x = SVector(0.0, 10_000.0, deg2rad(30), 0.0)
d = ThinDisc(0.0, Inf)

bins = collect(range(0.1, 1.5, 300))

eps = [-2.0, 0.0, 2.0]
data = map(eps) do ep
    m = JohannsenPsaltisMetric(M = 1.0, a = 0.4, ϵ3 = ep)
    g, f2 = lineprofile(m, x, d, verbose = true, bins = bins, maxrₑ=100.0)
end

data2 = map(eps) do ep
    m = JohannsenPsaltisMetric(M = 1.0, a = 0.6, ϵ3 = ep)
    g, f2 = lineprofile(m, x, d, verbose = true, bins = bins, maxrₑ=100.0)
end

begin
    fig = Figure(resolution = (480, 360))
    ga = fig[1,1] = GridLayout()
    ax = Axis(ga[2,1],
        ylabel = L"$F$ arb.",
        xlabel = L"E / E_\text{em}",
        title = "a = 0.4, θ = 30°"
    )
    ax2 = Axis(ga[2,2],
        xlabel = L"E / E_\text{em}",
        title = "a = 0.6, θ = 30°",
        yticks = [0.0, 0.5, 1.0],
    )
    hideydecorations!(ax2, grid=false)
    Makie.xlims!(ax, 0.28, 1.25)
    Makie.xlims!(ax2, 0.28, 1.25)

    norm = maximum(maximum(data[i][2]) for i in 1:3)
    lines = map(zip(data, eps)) do dat
        d, b = dat
        lines!(ax, d[1], d[2] ./ norm, label = Printf.@sprintf("b = %0.1f", b))
    end

    norm = maximum(maximum(data2[i][2]) for i in 1:3)
    for (dat) in zip(data2, eps)
        d, b = dat
        lines!(ax2, d[1], d[2] ./ norm, label = Printf.@sprintf("b = %0.1f", b))
    end

    leg = Legend(
        ga[1, 1:2],
        lines,
        ["ϵ = -2.0", "Kerr", "ϵ = 2.0"],
        orientation = :horizontal,
        # labelsize = 10,
        padding = (0, 0, 0, 0),
        framevisible = false,
    )
    rowgap!(ga, 5)
    colgap!(ga, 10)
    rowsize!(ga, 1, Auto(2))
    @savefigure(fig)
    fig
end 

