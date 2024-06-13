using CairoMakie, Makie, LaTeXStrings, Printf
using Gradus

include("common.jl")

function calculate_transfer_function(m, d, angle::Number, r0; offset = 1e-6, N = 500, kwargs...)
    x = @SVector [0.0, 1000.0, deg2rad(angle), 0.0]
    ctf = cunningham_transfer_function(m, x, d, r0, N = N; kwargs...)
    mask = @. (ctf.g✶ > offset) & (ctf.g✶ < 1 - offset)
    (;g = ctf.g✶[mask], f = ctf.f[mask], t = ctf.t[mask])
end

function tf_for_angle(m, d, angle; r = 4.0)
    high = @time calculate_transfer_function(m, d, angle, r)
    low = @time calculate_transfer_function(m, d, angle, r, N = 80)
    (;line = high, sample = low)
end

m = KerrMetric(M = 1.0, a = 0.998)
d = ThinDisc(0.0, 100.0)
tf1 = tf_for_angle(m, d, 74.0)
tf2 = tf_for_angle(m, d, 50.0)

m2 = KerrMetric(M = 1.0, a = 0.0)
tf3 = tf_for_angle(m2, d, 85.0; r = 11.0)
tf4 = tf_for_angle(m2, d, 50.0; r = 11.0)

begin
    palette = _default_palette()

    fig = Figure(resolution = (470, 330))
    ax1 = Axis(fig[2, 1], xlabel = L"g^\ast", ylabel = L"f", yticks = LinearTicks(4))

    c1 = popfirst!(palette)
    c2 = popfirst!(palette)
    c3 = popfirst!(palette)
    c4 = popfirst!(palette)

    lines!(ax1, tf1.line.g, tf1.line.f, color = c1, linewidth=1.5)
    scatter!(ax1, tf1.sample.g, tf1.sample.f, markersize = 7.0, color = c1)

    lines!(ax1, tf2.line.g, tf2.line.f, color = c2, linewidth=1.5)
    scatter!(ax1, tf2.sample.g, tf2.sample.f, markersize = 7.0, color = c2)

    lines!(ax1, tf3.line.g, tf3.line.f, color = c3, linewidth=1.5)
    scatter!(ax1, tf3.sample.g, tf3.sample.f, markersize = 7.0, color = c3)

    curves = [
      [LineElement(color = c, linestyle = nothing)
       MarkerElement(marker = :o, color = c)] for c in (c1, c2, c3)
    ]

    Legend(
        fig[1, 1],
        curves,
         [L"\theta=74^\circ", L"\theta=50^\circ",
         L"\theta=85^\circ"],
        orientation = :horizontal,
        framevisible = false,
        height = 10,
    )

    @savefigure(fig)
    fig
end
