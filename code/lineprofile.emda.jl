using Gradus
using Printf

include("common.jl")

x = SVector(0.0, 10_000.0, deg2rad(65), 0.0)
d = ThinDisc(0.0, Inf)

bins = collect(range(0.1, 1.5, 300))

bs = [0.0, 0.5, 1.0, 3.0]
data = map(bs) do b
    m = DilatonAxion(a = 0.7, b = b, β=0.0)
    g, f2 = lineprofile(m, x, d, verbose = true, bins = bins, maxrₑ=50.0)
end

x = SVector(0.0, 1000.0, deg2rad(25), 0.0)
data2 = map(bs) do b
    m = DilatonAxion(a = 0.7, b = b, β=0.0)
    g, f2 = lineprofile(m, x, d, verbose = true, bins = bins, maxrₑ=50.0)
end

begin
    fig = Figure(resolution = (480, 360))
    ga = fig[1,1] = GridLayout()
    ax = Axis(ga[2,1],
        ylabel = L"$F$ arb.",
        xlabel = L"E / E_\text{em}",
        title = "a = 0.7, θ = 65°"
    )
    ax2 = Axis(ga[2,2],
        xlabel = L"E / E_\text{em}",
        title = "a = 0.7, θ = 25°",
        yticks = [0.0, 0.5, 1.0],
    )
    hideydecorations!(ax2, grid=false)
    Makie.xlims!(ax, 0.3, 1.4)
    Makie.xlims!(ax2, 0.43, 1.13)

    norm = maximum(maximum(data[i][2]) for i in 1:4)
    lines = map(zip(data, bs)) do dat
        d, b = dat
        lines!(ax, d[1], d[2] ./ norm, label = Printf.@sprintf("b = %0.1f", b))
    end

    norm = maximum(maximum(data2[i][2]) for i in 1:4)
    for (dat) in zip(data2, bs)
        d, b = dat
        lines!(ax2, d[1], d[2] ./ norm, label = Printf.@sprintf("b = %0.1f", b))
    end

    leg = Legend(
        ga[1, 1:2],
        lines,
        ["Kerr", "b = 0.5", "b = 1.0", "b = 3.0"],
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
