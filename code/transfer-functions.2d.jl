using CairoMakie, Makie, LaTeXStrings, Printf
using Gradus

include("common.jl")

function transfer_branches(m, x, d, radii)
    d = ThinDisc(0.0, d.outer_radius)
    itb = @time Gradus.interpolated_transfer_branches(m, x, d, radii; verbose = true)
end

function calculate_2d_transfer_function(m, x, model, itb, prof, radii)
    bins = collect(range(0.0, 1.5, 500))
    tbins = collect(range(0, 100.0, 500))

    t0 = continuum_time(m, x, model)

    flux = @time Gradus.integrate_lagtransfer(
        prof,
        itb,
        radii,
        # range(4.0, 15.0, 2),
        bins,
        tbins;
        t0 = t0,
        Nr = 3000,
        h = 1e-8,
    )

    flux[flux.==0] .= NaN

    bins, tbins, flux
end


m = KerrMetric(M = 1.0, a = 0.998)
d = ThinDisc(Gradus.isco(m), 10000.0)
model = LampPostModel(h = 10.0, θ = deg2rad(0.01))

prof = @time emissivity_profile(m, d, model; n_samples = 1000)
radii = Gradus.Grids._inverse_grid(Gradus.isco(m), 150.0, 73)

x1 = SVector(0.0, 1e4, deg2rad(45), 0.0)
x2 = SVector(0.0, 1e4, deg2rad(80), 0.0)

itb1 = transfer_branches(m, x1, d, radii)
itb2 = transfer_branches(m, x2, d, radii)

E1, t1, f1 = calculate_2d_transfer_function(m, x1, model, itb1, prof, radii)
E2, t2, f2 = calculate_2d_transfer_function(m, x2, model, itb2, prof, radii)

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

    cmap = :batlow
    hm = heatmap!(ax1, t2, E2, log.(abs.(f2')), colormap = cmap, levels = 7)
    heatmap!(ax2, t1, E1, log.(abs.(f1')), colormap = cmap)

    # Colorbar(ga[1,2], hm)
    hidexdecorations!(ax1, grid = false)
    rowgap!(ga, 10)
    linkxaxes!(ax2, ax1)

    @savefigure(fig)
    fig
end
