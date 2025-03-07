using Gradus, Makie, CairoMakie, LaTeXStrings, Printf

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
    Gradus._threaded_map(eachrow(f)) do ψ
        t_extended, ψ_extended = Gradus.extend_domain_with_zeros(t, ψ, 1 / flo)
        Gradus.lag_frequency(t_extended, ψ_extended; kwargs...)
    end
end

function calculate_2d_transfer_function(m, x, model, itb, prof, radii)
    bins = collect(range(0.0, 1.7, 2000))
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
        g_grid_upscale = 1,
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

function calculate_for_disc(m, d)
    model = LampPostModel(h = 10.0)
    # different angles
    angles = [60, 80, 85]
    datas = map(angles) do angle
        @info "Angle = $angle"
        x = SVector(0.0, 10_000.0, deg2rad(angle), 0.0)
        b0 = angle < 45 ? 1.0 : 2.0
        itb = Gradus.interpolated_transfer_branches(m, x, d, radii; verbose = true, β₀ = b0, abstol = 1e-15, reltol = 1e-15)
        calculate_lag_transfer(m, x, d, model, radii, itb)
    end
    datas
end

m = KerrMetric(1.0, 0.998)
radii = Gradus.Grids._geometric_grid(Gradus.isco(m), 1000.0, 400)

begin
    @info "Thin disc"
    datas1 = calculate_for_disc(m, ThinDisc(0.0, Inf))
    @info "Thick disc 1"
    datas2 = calculate_for_disc(m, ShakuraSunyaev(m, eddington_ratio = 0.1))
    @info "Thick disc 2"
    datas3 = calculate_for_disc(m, ShakuraSunyaev(m, eddington_ratio = 0.2))
    @info "Thick disc 3"
    datas4 = calculate_for_disc(m, ShakuraSunyaev(m, eddington_ratio = 0.3))
end

begin
    palette = _default_palette()

    fig = Figure(size = (1100, 500))
    ga = fig[1, 1] = GridLayout()

    # frequency ranges for the energy calculations
    lims1 = (1e-3, 2e-3)
    lims2 = (4e-3, 8e-3)
    lims3 = (1.9e-2, 4e-2)
    lims4 = (1e-5, 1e-4)

    # x_lims1 = (0.83, 1.08)
    x_lims2 = (0.45, 1.45)
    x_lims1 = (0.45, 1.6)

    axes = map(enumerate((datas2, datas3, datas4))) do (col, ds)
        map(enumerate(ds)) do (row, d)
            ax = Axis(ga[row, col],
                ylabel ="τ",
                xlabel ="E/E₀",
                yticks = [0, 20, 40]
            )
            thin = datas1[row]

            for (lim) in ((lims1, lims2, lims3, lims4))
                thin_le = lag_energy(thin.all_pulses, lim...)
                lines!(ax, thin.E, thin_le, color = :grey, alpha=0.3)

                le = lag_energy(d.all_pulses, lim...)
                lines!(ax, d.E, le)
            end
            if col != 1
                hideydecorations!(ax, grid=false)
            end
            if row != 3
                hidexdecorations!(ax, grid=false)
            end
            xlims!(ax, x_lims1...)
            ax
        end
    end
    axs = reshape(reduce(vcat, axes), (length(axes[1]), length(axes)))
    for i in 1:length(axes)
        linkxaxes!(axs[1, i], axs[2, i])
        # linkxaxes!(axs[3, i], axs[4, i])
    end
    for i in 1:length(axes[1])
        root = axs[i, 1]
        for a in axs[i, 2:end]
            linkyaxes!(root, a)
        end
    end

    rowgap!(ga, 10)
    colgap!(ga, 10)

    for (i, a) in enumerate([0.1, 0.2, 0.3])
        Label(
            ga[1, i, Top()],
            text = "Ṁ / Ṁₑ = $a",
            padding = (0, 0, 10, 0),
            fontsize = 13,
            font = :bold,
        )
    end
    for (i, h) in enumerate([60, 80, 85])
        Label(
            ga[i, 3, Right()],
            text = "θ = $(h)°",
            padding = (10, 0, 00, 0),
            fontsize = 13,
            font = :bold,
        )
    end

    @savefigure(fig)
    fig
end
