using Gradus, Makie, CairoMakie, LaTeXStrings, Printf
import Plots

include("common.jl")

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

function calculate_2d_transfer_function(m, x, model, itb, prof, radii)
    bins = collect(range(0.0, 1.5, 2000))
    tbins = collect(range(0, 2000.0, 8000))

    t0 = continuum_time(m, x, model)
    @show t0

    flux = @time Gradus.integrate_lagtransfer(
        prof,
        itb,
        bins,
        tbins;
        t0 = t0,
        n_radii = 8000,
        rmin = minimum(radii),
        rmax = maximum(radii),
        h = 1e-8,
        g_grid_upscale = 3,
    )

    flux[flux.==0] .= NaN
    bins, tbins, flux
end

function calculate_lag_transfer(m, x, d, model, radii, itb)
    prof = @time emissivity_profile(m, d, model; n_samples = 150_000)
    E, t, f = @time calculate_2d_transfer_function(m, x, model, itb, prof, radii)
    ψ = Gradus.sum_impulse_response(f)
    freq, τ = @time lag_frequency(t, f)

    fn = replace(f, NaN => 0)
    all_pulses = lag_frequency_rowwise(t, fn ./ maximum(sum(fn, dims = 2)))
    (; freq, τ, ψ, t, all_pulses, E, f)
end

m = KerrMetric(1.0, 0.998)
angles = collect(5:10:85)
radii = Gradus.Grids._geometric_grid(Gradus.isco(m), 1000.0, 300)

model = LampPostModel(h = 10.0)
d = ShakuraSunyaev(m)

# different angles
# itbs = map(angles) do angle
#     @info "Angle = $angle"
#     x = SVector(0.0, 10_000.0, deg2rad(angle), 0.0)
#     b0 = angle < 45 ? 1.0 : 2.0
#     itb = Gradus.interpolated_transfer_branches(m, x, d, radii; verbose = true, β₀ = b0)
# end ;

# datas = map(zip(angles, itbs)) do it
#     angle, itb = it
#     @info "Angle = $angle"
#     x = SVector(0.0, 10_000.0, deg2rad(angle), 0.0)
#     freq1, τ1, impulse1, time1 = calculate_lag_transfer(m, x, d, model, radii, itb)
# end ;

# different rates
rates = [0.05, 0.1, 0.2, 0.3]
itbs = map(rates) do rate
    x = SVector(0.0, 10_000.0, deg2rad(60), 0.0)
    d = ShakuraSunyaev(m; eddington_ratio = rate)
    itb = Gradus.interpolated_transfer_branches(m, x, d, radii; verbose = true, β₀ = 2.0, abstol = 1e-12, reltol = 1e-12)
end ;

datas = map(zip(rates, itbs)) do it
    rate, itb = it
    x = SVector(0.0, 10_000.0, deg2rad(60), 0.0)
    d = ShakuraSunyaev(m; eddington_ratio = rate)
    calculate_lag_transfer(m, x, d, model, radii, itb)
end ;

begin
    for (n,dd) in zip(rates, datas)
        p = Plots.heatmap(dd.t, dd.E, log10.(dd.f), xlims = (0, 2000), title="$n")
        display(p)
    end
end

begin
    lims1 = (1e-3, 2e-3)
    lims2 = (4e-3, 8e-3)
    lims3 = (1.9e-2, 4e-2)
    lims4 = (1e-5, 1e-4)

    lagens = map(datas) do d
        le1 = lag_energy(d.all_pulses, lims1...)
        le2 = lag_energy(d.all_pulses, lims2...)
        le3 = lag_energy(d.all_pulses, lims3...)
        le4 = lag_energy(d.all_pulses, lims4...)
        (le1, le2, le3, le4)
    end
end

begin
    fig = Figure()
    ax = Axis(fig[1,1])

    for le in lagens[[1]]
        # lines!(ax, datas[1].E, le[1])
        # lines!(ax, datas[1].E, le[2])
        # lines!(ax, datas[1].E, le[3])
        lines!(ax, datas[1].E, le[4])
    end
    fig
end

begin
    palette = _default_palette()

    fig = Figure(resolution = (500, 650))
    ga = fig[1, 1] = GridLayout()
    ax1 =
        Axis(ga[1, 1], yscale = log10, xlabel = L"Time $t$", ylabel = L"Impulse Response$$")
    ax2 = Axis(
        ga[2, 1],
        xscale = log10,
        xminorgridvisible = true,
        xminorticks = IntervalsBetween(10),
        ylabel = L"Lag $\tau$",
        xlabel = L"Frequency $f$",
    )

    for dat in datas
        freq1, τ1, impulse1, time1 = dat
        c = popfirst!(palette)
        lines!(ax1, time1, abs.(impulse1), color = c)
        lines!(ax2, freq1, τ1, color = c)
    end

    xlims!(ax2, 5e-5, 0.3)
    ylims!(ax2, -8, 30)

    ylims!(ax1, 1e-3, 0.1)
    xlims!(ax1, 0, 50)

    Label(
        ga[1, 1, Right()],
        text = "a",
        padding = (8, 0, -90, 0),
        fontsize = 18,
        font = :bold,
    )
    Label(
        ga[2, 1, Right()],
        text = "b",
        padding = (8, 0, -410, 0),
        fontsize = 18,
        font = :bold,
    )

    # @savefigure(fig)
    fig
end
