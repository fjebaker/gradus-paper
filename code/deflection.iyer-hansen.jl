using Makie, CairoMakie, LaTeXStrings
using Gradus
using IyerHansen2009

function shadow(M, a; pos = true)
    s = (pos ? 1 : -1)
    -a + s * 6M * cos((1 / 3) * acos(-s * a / M))
end
shadow(m::AbstractMetric; kwargs...) = shadow(m.M, m.a; kwargs...)

m = KerrMetric(1.0, 0.9)
x = SVector(0.0, 100_000_000.0, π / 2, 0.0)

as = range(-20.0, 20.0, 300)
vs = [map_impact_parameters(m, x, a, 0.0) for a in as]
xs = fill(x, size(vs))

sols = tracegeodesics(
    m,
    xs,
    vs,
    300_000_000.0,
    save_on = false,
    chart = Gradus.chart_for_metric(m, 300_000_000.0),
    abstol = 1e-14,
    reltol = 1e-14,
)
points = unpack_solution(sols)

deviation = [
    if i.status != StatusCodes.WithinInnerBoundary
        i.x[4]
    else
        NaN
    end for i in points
]

calc = (abs.(deviation) .- π) ./ π

res = map(as) do a
    try
        deflection_angle(m.M, m.a, -a) ./ π
    catch
        NaN
    end
end


begin
    fig = Figure(resolution = (500, 420))
    ga = fig[1, 1] = GridLayout()
    ax = Axis(ga[1, 1], ylabel = "Deflection [rad]")
    ylims!(ax, 0.0, 2.0)

    poly!(
        ax,
        Point2f[
            (-shadow(m; pos = false), 0),
            (-shadow(m; pos = false), 2.0),
            (-shadow(m; pos = true), 2.0),
            (-shadow(m; pos = true), 0.0),
        ],
        color = RGBAf(0.1, 0.2, 0.3, 0.1),
        strokewidth = 1.0,
    )

    scatter!(ax, as, calc, marker = :x, color = :black, markersize = 11, label = L"\delta")
    lines!(ax, as, res, color = :red, label = L"\hat{\alpha}")
    axislegend(ax, position = :lt)

    ax2 = Axis(
        ga[2, 1],
        yscale = log10,
        yticks = LogTicks(LinearTicks(3)),
        xlabel = L"\alpha",
        ylabel = L"\text{abs}(\delta - \hat{\alpha})",
    )


    ax2_ymin = 2e-8
    ax2_ymax = 2e-5
    poly!(
        ax2,
        Point2f[
            (-shadow(m; pos = false), ax2_ymin),
            (-shadow(m; pos = false), ax2_ymax),
            (-shadow(m; pos = true), ax2_ymax),
            (-shadow(m; pos = true), ax2_ymin),
        ],
        color = RGBAf(0.1, 0.2, 0.3, 0.1),
        strokewidth = 1.0,
    )

    lines!(ax2, as, abs.(calc .- res), color = :black)
    ylims!(ax2, ax2_ymin, ax2_ymax)

    hidexdecorations!(ax, grid = false)
    linkxaxes!(ax2, ax)
    rowsize!(ga, 1, Auto(2))
    rowgap!(ga, 15)
    xlims!(ax2, -20, 20)

    fig
    @savefigure(fig)
end
