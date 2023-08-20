using CairoMakie, Makie, LaTeXStrings, Printf
using Gradus

include("common.jl")

function continuum_time(m, x, model)
    pos, _ = Gradus.sample_position_velocity(m, model)
    target = SVector(pos[2:end]...)
    _, _, gp, _ =
        Gradus.optimize_for_target(target, m, x, chart = Gradus.chart_for_metric(m, 2x[2]))
    @show gp.x[1]
end

function calculate_2d_transfer_function(m, x, d, model)
    prof = @time emissivity_profile(m, d, model; n_samples = 1000)
    radii = Gradus.Grids._inverse_grid(Gradus.isco(m), 150.0, 73)
    d = GeometricThinDisc(0.0, d.outer_radius, π / 2)
    itb = @time Gradus.interpolated_transfer_branches(m, x, d, radii; verbose = true)

    bins = collect(range(0.0, 1.5, 500))
    tbins = collect(range(0, 100.0, 500))

    t0 = continuum_time(m, x, model)

    flux = @time Gradus.integrate_lagtransfer(
        prof,
        itb,
        radii,
        bins,
        tbins;
        t0 = t0,
        Nr = 3000,
        h = 1e-12,
    )

    flux[flux.==0] .= NaN

    bins, tbins, flux
end


m = KerrMetric(M = 1.0, a = 0.998)
d = GeometricThinDisc(Gradus.isco(m), 10000.0, π / 2)
model = LampPostModel(h = 10.0, θ = deg2rad(0.01))

x = SVector(0.0, 1e6, deg2rad(45), 0.0)
E1, t1, f1 = calculate_2d_transfer_function(m, x, d, model)
x = SVector(0.0, 1e6, deg2rad(80), 0.0)
E2, t2, f2 = calculate_2d_transfer_function(m, x, d, model)

begin
    fig = Figure(resolution = (500, 550))
    ga = fig[1, 1] = GridLayout()

    ax1 = Axis(
        ga[1, 1],
        ylabel = L"E / E_\text{line}",
        xlabel = L"t",
        xticks = LinearTicks(5),
    )
    ylims!(ax1, nothing, 1.52)
    ax2 = Axis(
        ga[2, 1],
        ylabel = L"E / E_\text{line}",
        xlabel = L"t",
        xticks = LinearTicks(5),
    )
    ylims!(ax2, nothing, 1.2)

    hm = heatmap!(ax1, t2, E2, log.(abs.(f2')))
    heatmap!(ax2, t1, E1, log.(abs.(f1')))

    # Colorbar(ga[1,2], hm)
    hidexdecorations!(ax1, grid = false)
    rowgap!(ga, 10)
    linkxaxes!(ax2, ax1)

    fig
    @savefigure(fig)
end
