using Makie, CairoMakie, LaTeXStrings
using Gradus

include("common.jl")

function trace_even(m, d, model)
    δs = deg2rad.(range(0.01, 179.99, 20))

    x, v = Gradus.sample_position_velocity(m, model)
    vs = map(δs) do δ
        Gradus.sky_angles_to_velocity(m, x, v, δ, π)
    end
    xs = fill(x, size(vs))

    sols = tracegeodesics(m, xs, vs, d, 2000.0)
end

function plot_paths_xz!(ax, sols; N = 700, kwargs...)
    for sol in sols
        x, y, z = Gradus._extract_path(sol, N)
        lines!(ax, x, z; kwargs...)
    end
end

m = KerrMetric(a = 0.998)
d = GeometricThinDisc(Gradus.isco(m), 1000.0, π / 2)

model = LampPostModel(θ = 1e-6)
sols = tracegeodesics(
    m,
    model,
    d,
    2000.0,
    n_samples = 64,
    sampler = EvenSampler(BothHemispheres(), GoldenSpiralGenerator()),
)

sols1 = trace_even(m, d, model)

sols_all = vcat(sols1.u[1:8], sols1.u[10:end]);
sols_selected = sols1.u[8:9]


R = Gradus.inner_radius(m)
angles = collect(range(0.0, π, 100))

horizon_x = @. cos(angles) * R
horizon_y = @. sin(angles) * R

begin
    fig = Figure(resolution = (500, 500))
    ga = fig[1, 1] = GridLayout()

    ax = Axis(
        ga[1, 1],
        aspect = DataAspect(),
        yticks = LinearTicks(3),
        ylabel = L"z\, [\rg]",
    )
    xlims!(ax, -10, 10)
    ylims!(ax, 0, 10)
    lines!(ax, horizon_x, horizon_y, color = :black, linewidth = 3.0)
    plot_paths_xz!(ax, sols)


    ax2 = Axis(
        ga[2, 1],
        aspect = DataAspect(),
        yticks = LinearTicks(3),
        xlabel = L"x\,[\rg]",
        ylabel = L"z\, [\rg]",
    )
    xlims!(ax2, -10, 10)
    ylims!(ax2, 0, 10)
    lines!(ax2, horizon_x, horizon_y, color = :black, linewidth = 3.0)
    plot_paths_xz!(ax2, sols1)

    linkxaxes!(ax, ax2)
    hidexdecorations!(ax, grid = false)
    rowgap!(ga, 10)

    ax3 = Axis(
        ga[2, 1],
        width = Relative(0.35),
        height = Relative(0.5),
        halign = 0.1,
        aspect = DataAspect(),
    )
    xlims!(ax3, 0.00, 0.002)
    ylims!(ax3, 4.9992, 5.0)
    hidedecorations!(ax3)

    plot_paths_xz!(ax3, sols_all, color = :lightgray)
    plot_paths_xz!(ax3, sols_selected)

    text!(ax3, (0.0016, 4.99942), text = L"\Delta \theta")
    text!(ax2, sols_selected[1].u[end][2], 0.0, text = L"\Delta r")

    Label(ga[1, 1, Right()], "a", padding = (10, 00, 170, 0), font = :bold, fontsize = 20)
    Label(ga[2, 1, Right()], "b", padding = (10, 00, 170, 0), font = :bold, fontsize = 20)

    fig
    #@savefigure(fig)
end
