using CairoMakie, Makie, LaTeXStrings
using Printf
using Gradus

include("common.jl")

function E_Lz(m::AbstractMetric, r; q = 0.0)
    CircularOrbits.energy_angmom(m, SVector(r, π / 2))
end

function el_for_metric(m, r_targets)
    tol = 1e-12
    vels = Gradus.solve_equatorial_circular_orbit(
        m,
        r_targets,
        max_time = 30.0,
        abstol = tol,
        reltol = tol,
        lower_rate =  0.99
    )
    
    stabilities = map(
        i -> Gradus.measure_stability(m, r_targets[i], vels[i]),
        eachindex(r_targets),
    )
    @show extrema(stabilities)

    t_final = 2000.0
    sols = map(zip(r_targets, vels)) do (r, v)
        Gradus.trace_single_orbit(
            m,
            r,
            v,
            max_time = t_final,
            abstol = tol,
            reltol = tol,
        )
    end
    map(sols) do sol
        v = sol.u[1][5:8]
        x = sol.u[1][1:4]

        rtheta = x[2:3]
        g = Gradus.metric(m, rtheta)
        utuphi = (g*v)[[1, 4]]
        E = CircularOrbits.energy(m, rtheta, utuphi)
        L = CircularOrbits.angmom(m, rtheta, utuphi)

        (E, L)
    end
end

function plot_curve!(ax, m, rs; q = 0.0, kwargs...)
    values = E_Lz.(m, rs; q)
    Es = first.(values)
    Lzs = last.(values)
    lines!(ax, Lzs, Es; kwargs...)
end

logspace(start, stop, num::Int) = collect(Gradus.Grids._geometric_grid(start, stop, num))

m0 = KerrMetric(1.0, 0.0)
rs0 = collect(range(3.3, 50.0, 300))
comp0 = el_for_metric(m0, logspace(3.3, 50.0, 30))

m1 = KerrMetric(1.0, 0.4)
rs1 = collect(range(2.5, 50.0, 300))
comp1 = el_for_metric(m1, logspace(2.5, 50.0, 30))

m2 = KerrMetric(1.0, 0.8)
rs2 = collect(range(1.96, 50.0, 300))
comp2 = el_for_metric(m2, logspace(1.99, 50.0, 30))

m3 = KerrMetric(1.0, 0.998)
rs3 = collect(Gradus.Grids._geometric_grid(1.08, 50.0, 300))
comp3 = el_for_metric(m3, logspace(1.4, 50.0, 30))

begin
    fig = Figure(resolution = (450, 300))
    layout = fig[1, 1] = GridLayout()
    ax = Axis(layout[2, 1], xlabel = L"L_z", ylabel = L"E", yticks = LinearTicks(3))
    ylims!(ax, 0.65, 1.15)
    xlims!(ax, 1, 6.5)

    palette = _default_palette()
    s0 = plot_curve!(ax, m0, rs0, color = popfirst!(palette))
    k1 = plot_curve!(ax, m1, rs1, linewidth = 1.5, color = popfirst!(palette))
    k2 = plot_curve!(ax, m2, rs2, linewidth = 1.5, color = popfirst!(palette))
    k3 = plot_curve!(ax, m3, rs3, linewidth = 1.5, color = popfirst!(palette))
    
    palette = _default_palette()
    scatter_kwargs = (; markersize = 8, marker = :x)
    scatter!(ax, last.(comp0), first.(comp0); color = popfirst!(palette), scatter_kwargs...)
    scatter!(ax, last.(comp1), first.(comp1); color = popfirst!(palette), scatter_kwargs...)
    scatter!(ax, last.(comp2), first.(comp2); color = popfirst!(palette), scatter_kwargs...)
    scatter!(ax, last.(comp3), first.(comp3); color = popfirst!(palette), scatter_kwargs...)

    a1_str = Printf.@sprintf("%0.1f", m1.a)
    a2_str = Printf.@sprintf("%0.1f", m2.a)
    a3_str = Printf.@sprintf("%0.1f", m3.a)

    ga = layout[1, 1] = GridLayout()
    Legend(
        ga[1, 1],
        [k3, k2, k1, s0],
        [L"a = %$(a)" for a in reverse(("0.0", a1_str, a2_str, a3_str))],
        orientation = :horizontal,
        framevisible = false,
        height = 10,
    )
    colgap!(ga, 0)
    rowgap!(ga, 0)

    rowgap!(layout, 5)
    colgap!(layout, 0)

    resize_to_layout!(fig)
    @savefigure(fig)
    fig
end
