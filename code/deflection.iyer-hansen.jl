using Makie, CairoMakie, LaTeXStrings
using Gradus
using IyerHansen2009

function trace_deflection(m::AbstractMetric, xs, vs; kwargs...)
    sols = tracegeodesics(
        m,
        xs,
        vs,
        400_000_000.0,
        save_on = false,
        chart = Gradus.chart_for_metric(m, 400_000_000.0),
        abstol = 1e-14,
        reltol = 1e-14;
        kwargs...,
    )
    points = unpack_solution(sols)
    deflection = [
        if i.status != StatusCodes.WithinInnerBoundary
            i.x[4]
        else
            NaN
        end for i in points
    ]
    (abs.(deflection) .- π)
end

function anayltic_deflection(m::AbstractMetric, a)
    try
        deflection_angle(m.M, m.a, -a)
    catch
        NaN
    end
end

function shadow(M, a; pos = true)
    s = (pos ? 1 : -1)
    -a + s * 6M * cos((1 / 3) * acos(-s * a / M))
end
shadow(m::AbstractMetric; kwargs...) = shadow(m.M, m.a; kwargs...)

m = KerrMetric(1.0, 0.9)
x = SVector(0.0, 200_000_000.0, π / 2, 0.0)

as = range(-20.0, 20.0, 300)
vs = [map_impact_parameters(m, x, a, 0.0) for a in as]
xs = fill(x, size(vs))

analytic = map(i -> anayltic_deflection(m, i), as)

calc = trace_deflection(m, xs, vs)
calc_2 = trace_deflection(m, xs, vs; solver = Gradus.Feagin10())
calc_3 = trace_deflection(m, xs, vs; solver = Gradus.Vern6())
calc_4 = trace_deflection(m, xs, vs; solver = Gradus.RK4())

begin
    fig = Figure(resolution = (500, 460))
    ga = fig[1, 1] = GridLayout()
    ax = Axis(ga[1, 1], ylabel = L"|Deflection| $/\pi$")
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

    scatter!(ax, as, calc ./ π, markersize = 6, label = L"\delta")
    lines!(ax, as, analytic ./ π, color = :black, label = L"\hat{\alpha}")
    axislegend(ax, position = :lt)

    ax2 = Axis(
        ga[3, 1],
        yscale = log10,
        yticks = LogTicks(LinearTicks(3)),
        xlabel = L"\alpha",
        ylabel = L"\text{abs}(\delta x^\phi - \hat{\alpha})",
    )


    ax2_ymin = 2e-10
    ax2_ymax = 2e-2
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

    l1 = lines!(ax2, as, abs.(calc .- analytic))
    l2 = lines!(ax2, as, abs.(calc_2 .- analytic))
    l3 = lines!(ax2, as, abs.(calc_3 .- analytic))
    l4 = lines!(ax2, as, abs.(calc_4 .- analytic))
    ylims!(ax2, ax2_ymin, ax2_ymax)

    leg = Legend(
        ga[2, 1],
        [l1, l2, l3, l4],
        ["Tsit5", "Feagin10", "Vern6", "RK4"],
        orientation = :horizontal,
        # labelsize = 10,
        padding = (0, 0, 0, 0),
        framevisible = false,
    )

    hidexdecorations!(ax, grid = false)
    linkxaxes!(ax2, ax)
    rowsize!(ga, 1, Auto(2))
    rowgap!(ga, 15)
    xlims!(ax2, -20, 20)

    resize_to_layout!(fig)
    fig
    # @savefigure(fig)
end
