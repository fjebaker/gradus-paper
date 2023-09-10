using CairoMakie, Makie, LaTeXStrings
using Gradus
using Relxill
using SpectralFitting

include("common.jl")

function evaluate_relline(energy, m::AbstractMetric, x, maxrₑ)
    flux_relline = invokemodel(
        energy,
        XS_Relline(
            lineE = FitParam(1.0),
            a = FitParam(m.a),
            outer_r = FitParam(maxrₑ),
            θ_obs = FitParam(rad2deg(x[3])),
            limb = FitParam(0.0),
        ),
    )
    flux_relline
end

# emissivity function
ε(r) = r^(-3)

function compute_lineprofiles(gs, m::KerrMetric, x)
    d = ThinDisc(0.0, 500.0)
    redshift = ConstPointFunctions.redshift(m, x)
    # maximal integration radius
    maxrₑ = 50.0

    # g grid to do flux integration over
    _, flux = @time lineprofile(
        gs,
        ε,
        m,
        x,
        d,
        redshift_pf = redshift,
        maxrₑ = maxrₑ,
        verbose = true,
    )
    # transform to observed energy
    flux_relline = evaluate_relline(gs, m, x, maxrₑ)

    flux, flux_relline
end

x = SVector(0.0, 1000.0, deg2rad(40), 0.0)

m = KerrMetric(M = 1.0, a = 0.998)
g1 = range(0.05, 1.2, 500) |> collect
flux, relline_flux = compute_lineprofiles(g1, m, x)

m = KerrMetric(1.0, 0.0)
g2 = range(0.4, 1.2, 500) |> collect
flux2, relline_flux2 = compute_lineprofiles(g2, m, x)

begin
    norm_factor = maximum(flux)
    f1_g = flux ./ norm_factor
    f1_r = relline_flux ./ norm_factor
end
begin
    norm_factor = maximum(flux2)
    f2_g = flux2 ./ norm_factor
    f2_r = relline_flux2 ./ norm_factor
end

function plot_lineprofile(ax1, ax2, x, f1, f2, c1, c2)
    lines!(ax1, x, f1, label = "Gradus.jl", color = c1)
    lines!(ax1, x[1:end-1], f2, label = "relline", color = c2)

    lines!(ax2, x[1:end-1], f1[1:end-1] .- f2)

    hidexdecorations!(ax1, grid = false)

    ylims!(ax2, -0.015, 0.015)
    xlims!(ax2, extrema(x)...)
    linkxaxes!(ax2, ax1)
end

begin
    palette = Iterators.Stateful(Iterators.Cycle(Makie.wong_colors()))
    fig = Figure(resolution = (530, 400))

    ga = fig[1, 1] = GridLayout()

    color_gradus = popfirst!(palette)
    color_relline = popfirst!(palette)

    axf1 = Axis(ga[1, 1], ylabel = L"$F$ arb.", title = L"a = 0.998")
    axd1 = Axis(
        ga[2, 1],
        yticks = LinearTicks(3),
        xlabel = L"E / E_\text{em}",
        ylabel = L"\Delta F",
    )

    axf2 = Axis(ga[1, 2], title = L"a=0")
    axd2 = Axis(ga[2, 2], yticks = LinearTicks(3), xlabel = L"E / E_\text{em}")

    hideydecorations!(axf2, grid = false)
    hideydecorations!(axd2, grid = false)
    linkyaxes!(axd1, axd2)
    linkyaxes!(axf1, axf2)

    plot_lineprofile(axf1, axd1, g1, f1_g, f1_r, color_gradus, color_relline)
    plot_lineprofile(axf2, axd2, g2, f2_g, f2_r, color_gradus, color_relline)

    axislegend(axf1, position = (:left, :top), framevisible = false)

    rowsize!(ga, 1, Auto(3))
    colgap!(ga, 10)
    rowgap!(ga, 15)

    resize_to_layout!(fig)
    @savefigure(fig)
    fig
end
