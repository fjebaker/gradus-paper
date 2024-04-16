using Makie, CairoMakie, LaTeXStrings
using Gradus

include("common.jl")

function calculate_2d_transfer_function(m, x, model, itb, prof, radii)
    bins = collect(range(0.0, 1.5, 1000))
    tbins = collect(range(0, 2000.0, 3000))

    t0 = continuum_time(m, x, model)
    @show t0

    flux = @time Gradus.integrate_lagtransfer(
        prof,
        itb,
        radii,
        bins,
        tbins;
        t0 = t0,
        Nr = 8000,
        h = 1e-8,
        g_grid_upscale = 3,
    )

    flux[flux.==0] .= NaN

    bins, tbins, flux
end

function lag_energy(data, flow, fhi)
    f_example = data[1][1]
    i1 = findfirst(>(flow), f_example)
    i2 = findfirst(>(fhi), f_example)
    lE = map(data) do (_, tau)
        sum(tau[i1:i2]) / (i2 - i1)
    end
end

function lag_frequency_rowwise(t, f::AbstractMatrix; flo = 1e-5, kwargs...)
    map(eachrow(f)) do ψ
        t_extended, ψ_extended = Gradus.extend_domain_with_zeros(t, ψ, 1 / flo)
        Gradus.lag_frequency(t_extended, ψ_extended; kwargs...)
    end
end


m = KerrMetric(1.0, 0.998)
x = SVector(0.0, 10_000.0, deg2rad(45), 0.0)
radii = Gradus.Grids._inverse_grid(Gradus.isco(m), 1000.0, 300)

model = LampPostModel(h = 10.0)

# thin disc
d = ThinDisc(0.0, Inf)
itb = Gradus.interpolated_transfer_branches(m, x, d, radii; verbose = true, β₀ = 2.0)
prof = @time emissivity_profile(m, d, model; n_samples = 100_000)

E, t, f = @time calculate_2d_transfer_function(m, x, model, itb, prof, radii)
freq, τ = @time lag_frequency(t, f)

thick_d = ShakuraSunyaev(m)
thick_itb =
    Gradus.interpolated_transfer_branches(m, x, thick_d, radii; verbose = true, β₀ = 2.0)
thick_prof = @time emissivity_profile(m, thick_d, model; n_samples = 100_000)

thick_E, thick_t, thick_f =
    @time calculate_2d_transfer_function(m, x, model, thick_itb, thick_prof, radii)
thick_freq, thick_τ = @time lag_frequency(thick_t, thick_f)

begin
    data = lag_frequency_rowwise(t, replace(f, NaN => 0))
    thick_data = lag_frequency_rowwise(thick_t, replace(thick_f, NaN => 0))

    lims1 = (1e-3, 2e-3)
    lims2 = (4e-3, 8e-3)
    lims3 = (1.9e-2, 4e-2)
    lims4 = (1e-5, 1e-4)

    le1 = lag_energy(data, lims1...)
    le2 = lag_energy(data, lims2...)
    le3 = lag_energy(data, lims3...)
    le4 = lag_energy(data, lims4...)

    thick_le1 = lag_energy(thick_data, lims1...)
    thick_le2 = lag_energy(thick_data, lims2...)
    thick_le3 = lag_energy(thick_data, lims3...)
    thick_le4 = lag_energy(thick_data, lims4...)
end

begin
    palette = _default_palette()
    fig = Figure(resolution = (500, 400))
    ga = fig[1, 1] = GridLayout()
    ax = Axis(ga[1, 1], xlabel = L"E / E_\text{em}", ylabel = L"\tau")

    axmini = Axis(
        ga[1, 1],
        width = Relative(0.5),
        height = Relative(0.4),
        halign = 0.25,
        valign = 0.9,
        xscale = log10,
        xminorticks = IntervalsBetween(10),
        xminorgridvisible = true,
        backgroundcolor = :white,
        xlabel = L"f",
        ylabel = L"\tau",
    )
    translate!(axmini.scene, 0, 0, 10)
    # this needs separate translation as well, since it's drawn in the parent scene
    translate!(axmini.elements[:background], 0, 0, 9)
    translate!(axmini.elements[:xgridlines], 0, 0, 9)
    translate!(axmini.elements[:xminorgridlines], 0, 0, 9)
    translate!(axmini.elements[:ygridlines], 0, 0, 9)

    color1 = popfirst!(palette)
    color2 = popfirst!(palette)
    color3 = popfirst!(palette)
    color4 = popfirst!(palette)

    vspan!(
        axmini,
        [lims1[1], lims2[1], lims3[1], lims4[1]],
        [lims1[2], lims2[2], lims3[2], lims4[2]],
        color = [(c, 0.4) for c in [color1, color2, color3, color4]],
    )
    lines!(axmini, freq, τ, color = color1)
    lines!(axmini, thick_freq, thick_τ, linestyle = :dot, color = color1)

    lines!(ax, E, le1, color = color1)
    lines!(ax, E, le2, color = color2)
    lines!(ax, E, le3, color = color3)
    lines!(ax, E, le4, color = color4)

    lines!(ax, thick_E, thick_le1, linestyle = :dot, color = color1)
    lines!(ax, thick_E, thick_le2, linestyle = :dot, color = color2)
    lines!(ax, thick_E, thick_le3, linestyle = :dot, color = color3)
    lines!(ax, thick_E, thick_le4, linestyle = :dot, color = color4)


    xlims!(axmini, 5e-5, 0.3)
    xlims!(ax, 0.25, 1.23)

    @savefigure(fig)
    fig
end
