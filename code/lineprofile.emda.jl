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

begin
    fig = Figure(resolution = (420, 370))
    ax = Axis(fig[1,1],
        ylabel = L"$F$ arb.",
        xlabel = L"E / E_\text{em}",
        title = "a = 0.7, θ = 65°"
    )
    Makie.xlims!(ax, 0.3, 1.4)
    norm = maximum(data[1][2])
    lines = map(zip(data, bs)) do dat
        d, b = dat
        lines!(d[1], d[2] ./ norm, linewidth=2.0, label = Printf.@sprintf("b = %0.1f", b))
    end
    axislegend(ax, position = (:left, :top), framevisible = false)
    @savefigure(fig)
    fig
end 
