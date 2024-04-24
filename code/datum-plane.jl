using CairoMakie, Makie, LaTeXStrings, Printf
using Gradus

function traceimpact(m, x, α, β; kwargs...)
    v = map_impact_parameters(m, x, α, β)
    tracegeodesics(m, x, v, 2 * x[2]; kwargs...)
end

function traceimpact(m, x, α, β, d::Gradus.AbstractAccretionGeometry; kwargs...)
    v = map_impact_parameters(m, x, α, β)
    tracegeodesics(m, x, v, d, 2 * x[2]; kwargs...)
end

function plot_radius!(ax, r; kwargs...)
    ϕs = range(0.0, 2π, 100)
    xs = @. sin(ϕs) * r
    ys = @. cos(ϕs) * r
    lines!(ax, xs, ys; kwargs...)
end

include("common.jl")

m = KerrMetric(M = 1.0, a = 0.69)
d = Gradus.ShakuraSunyaev(m; eddington_ratio = 0.68)

plane = Gradus.DatumPlane(2.9)

rs = collect(range(Gradus.isco(m), 100.0, 300))
zs = Gradus.cross_section.(d, rs)

x = SVector(0.0, 1000.0, deg2rad(60), 0.0)

function draw_axis_key!(ax; x = 0.0, y = 0.0,  x_len = 2, y_len = 2, kwargs...)
    lines!(ax, [x, x, x + x_len], [y + y_len, y, y]; kwargs...)
    text!(ax, x - 0.4, y + y_len, text = L"z")
    text!(ax, x + x_len + 0.1, y - 0.5, text = L"r")
end

begin
    palette = _default_palette()

    fig = Figure(resolution = (470, 250))
    ax1 = Axis(
        fig[1, 1],
        xlabel = L"g^\ast",
        ylabel = L"f",
        yticks = LinearTicks(4),
        aspect = DataAspect(),
        xgridvisible = false,
        ygridvisible = false,
        xlabelvisible = false,
        ylabelvisible = false,
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0, 0),
        topspinevisible = false,
        leftspinevisible = false,
        rightspinevisible = false,
        bottomspinevisible = false,
    )
    hidedecorations!(ax1)

    r_em = 4.80
    vlines!(ax1, [-r_em, r_em], color = :black, linestyle = :dash, linewidth = 1.0)
    text!(ax1, r_em - 2, 15, text = L"r_\text{em}")
    text!(ax1, 20.0, 0, text = L"z=0")

    chosen_range = [7.0, 0.0, -2.0, 2.0]
    chosen_range = [7.0, 0.0]

    for b in range(-5, 10.0, step = 1)
        c = :lightgray
        sol = traceimpact(m, x, 0.0, b, d)
        X, _, Z = Gradus._extract_path(sol, 500, t_span = 60)
        lines!(ax1, X, Z, color = c, linewidth = 0.5)
    end

    plot_radius!(ax1, Gradus.inner_radius(m), color = :black, linewidth = 3.0)
    lines!(ax1, rs, zs, color = :black)
    lines!(ax1, rs, -zs, color = :black)
    lines!(ax1, -rs, -zs, color = :black)
    lines!(ax1, -rs, zs, color = :black)

    palette = _default_palette()
    hlines!(ax1, [plane.height], linewidth = 2.0, color = popfirst!(palette))
    hlines!(ax1, [0], linewidth = 1.0, color = :black, linestyle = :dash)

    for b in chosen_range
        c = popfirst!(palette)
        sol = traceimpact(m, x, 0.0, b, d)
        X, _, Z = Gradus._extract_path(sol, 500, t_span = 60)
        lines!(ax1, X, Z, color = c)
        scatter!(ax1, [X[end]], [Z[end]], color = c)
        sol = traceimpact(m, x, 0.0, b, plane)
        X, _, Z = Gradus._extract_path(sol, 500, t_span = 60)
        lines!(ax1, X, Z, linestyle = :dash, color = c)
        scatter!(ax1, [X[end]], [Z[end]], color = c)
    end

    draw_axis_key!(ax1, x = -9, y =14, color = :black)

    Label(
      fig[1,1,Top()],
        text = "A",
        padding = (125, 0, 0, 0),
        fontsize = 13,
        font = :bold,
    )

    Label(
      fig[1,1,Top()],
        text = "B",
        padding = (420, 0, 0, 0),
        fontsize = 13,
        font = :bold,
    )

    x_lim = 20
    offset = 10.0
    xlims!(ax1, -x_lim + offset, x_lim + offset)
    ylims!(ax1, -2.0, 17)

    @savefigure(fig)
    fig
end
