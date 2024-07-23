using CairoMakie, Makie, LaTeXStrings
using Gradus

include("common.jl")

function target_trajectory(m, x, pos)
    target = SVector(pos[2:end]...)
    _, _, gp, accuracy = Gradus.optimize_for_target(
        target,
        m,
        x;
        chart = Gradus.chart_for_metric(m, 2x[2]),
        callback = domain_upper_hemisphere(),
    )
    @show pos[2:3], gp.x[2:3], accuracy
    tracegeodesics(m, x, gp.v_init, ThinDisc(0.0, Inf), gp.λ_max)
end

function target_trajectory(m, x, model::Gradus.AbstractCoronaModel)
    pos, _ = Gradus.sample_position_velocity(m, model)
    target_trajectory(m, x, pos)
end

xfm(x) = (x)

function plot_circle!(ax, r; kwargs...)
    theta = collect(range(0, 2π, 200))
    x = @. xfm(r) * cos(theta)
    y = @. xfm(r) * sin(theta)
    lines!(ax, x, y; kwargs...)
end

function trace(m, x, model)
    NN = 256
    d = ThinDisc(0.0, Inf)
    sol1 = target_trajectory(m, x, model)
    coords1 = Gradus._extract_path(sol1, NN, t_span = 100.0, projection = :polar)

    n_samples = 1024
    begin
        δs = -deg2rad.(range(0.01, 179.99, n_samples))
        # we assume a point source
        xsrc, vsrc = Gradus.sample_position_velocity(m, model)
        velfunc = Gradus.polar_angle_to_velfunc(m, xsrc, vsrc, δs)
        gps = tracegeodesics(
            m,
            xsrc,
            velfunc,
            d,
            2000.0,
            save_on = false,
            ensemble = Gradus.EnsembleEndpointThreads(),
            trajectories = length(δs),
        )
        _, i = findmin(i -> i.x[1], gps)
        gp_min = gps[i]
    end

    sol2 = tracegeodesics(m, gp_min.x_init, gp_min.v_init, d, gp_min.λ_max)
    coords2 = Gradus._extract_path(sol2, NN, t_span = 100.0, projection = :polar)

    sol3 = target_trajectory(m, x, gp_min.x)
    coords3 = Gradus._extract_path(sol3, NN, t_span = 100.0, projection = :polar)

    logXYZ(coords1...), logXYZ(coords2...), logXYZ(coords3...)
end

function logXYZ(r, theta, phi)
    R = xfm.(r)
    x = @. R * sin(theta) * cos(phi)
    y = @. R * sin(theta) * sin(phi)
    z = @. R * cos(theta)
    (x, y, z)
end

m = KerrMetric(1.0, 0.95)
x = SVector(0.0, 29.0, deg2rad(45), 0.0)

model1 = LampPostModel(h = 3.3, θ = 0.001)
model2 = LampPostModel(h = 11.0, θ = 0.001)

(x1, y1, z1), (dx1, dy1, dz1), (dox1, doy1, doz1) = trace(m, x, model1)
(x2, y2, z2), (dx2, dy2, dz2), (dox2, doy2, doz2) = trace(m, x, model2)

begin
    palette = _default_palette()
    fig = Figure(resolution = (380, 320))
    ax = Axis(
        fig[1, 1],
        aspect = DataAspect(),
        xgridvisible = false,
        ygridvisible = false,
        xlabelvisible = false,
        ylabelvisible = false,
        leftspinevisible = false,
        topspinevisible = false,
        bottomspinevisible = false,
        rightspinevisible = false,
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0, 0),
    )
    Makie.hidedecorations!(ax)

    c2 = popfirst!(palette)
    c1 = popfirst!(palette)

    _ = popfirst!(palette)
    plot_circle!(ax, Gradus.inner_radius(m), color = :black, linewidth = 3.0)
    R = 18
    plot_circle!(ax, R, linestyle = :dash, color = popfirst!(palette))

    lines!(ax, x1, z1, color = c1, linestyle = :dot)
    lines!(ax, dx1, dz1, color = c1)
    lines!(ax, dox1, doz1, color = c1)

    lines!(ax, x2, z2, color = c2, linestyle = :dot)
    lines!(ax, dx2, dz2, color = c2)
    lines!(ax, dox2, doz2, color = c2)

    scatter!(ax, [last(x1)], [last(z1)], color = c1)
    scatter!(ax, [last(x2)], [last(z2)], color = c2)

    xlims!(ax, -3, 25)
    ylims!(ax, -2, 23)

    text!(ax, [last(x1) - 1.8], [last(z1) + 0.7], text = L"h_1", fontsize = 18)
    text!(ax, [last(x2) - 1.8], [last(z2) + 0.7], text = L"h_2", fontsize = 18)

    theta = deg2rad(15)
    text!(ax, [R * cos(theta) + 0.6], [R * sin(theta)], text = L"R", fontsize = 18)

    text!(
        ax,
        [x[2] * sin(x[3]) - 3.5],
        [x[2] * cos(x[3])],
        text = L"x_\text{obs}",
        fontsize = 18,
    )
    draw_observer_eye!(
        ax,
        [x[2] * sin(x[3]) + 1.6],
        [x[2] * cos(x[3]) + 0.6 + 1],
        1,
        flip = true,
        linewidth = 1.7,
        rot = -x[3],
    )

    lines!(ax, [Gradus.isco(m), 100.0], [0.0, 0.0], color = :black, linewidth = 2.0)
    lines!(ax, [-Gradus.isco(m), -100.0], [0.0, 0.0], color = :black, linewidth = 2.0)

    @savefigure(fig)
    fig
end
