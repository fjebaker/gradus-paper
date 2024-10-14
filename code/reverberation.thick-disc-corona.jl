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
        g_grid_upscale = 3,
    )

    flux[flux.==0] .= NaN
    bins, tbins, flux
end

function calculate_lag_transfer(m, x, d, model, radii, itb)
    prof = @time emissivity_profile(m, d, model; n_samples = 200_000)
    E, t, f = @time calculate_2d_transfer_function(m, x, model, itb, prof, radii)
    ψ = Gradus.sum_impulse_response(f)
    freq, τ = @time lag_frequency(t, f)

    fn = replace(f, NaN => 0)
    all_pulses = lag_frequency_rowwise(t, fn ./ maximum(sum(fn, dims = 2)))
    (; freq, τ, ψ, t, all_pulses, E, f)
end

function calculate_for_model(m, d, model)
    # different angles
    angles = [10, 30, 60, 80]
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
    @info "Thick disc 1"
    m_datas2 = calculate_for_model(m, ShakuraSunyaev(m, eddington_ratio = 0.3), LampPostModel(h=2.0))
    @info "Thick disc 2"
    m_datas3 = calculate_for_model(m, ShakuraSunyaev(m, eddington_ratio = 0.3), LampPostModel(h=5.0))
    @info "Thick disc 3"
    m_datas4 = calculate_for_model(m, ShakuraSunyaev(m, eddington_ratio = 0.3), LampPostModel(h=10.0))
end

begin
    @info "Thin disc 1"
    t_datas2 = calculate_for_model(m, ThinDisc(0.0, Inf), LampPostModel(h=2.0))
    @info "Thin disc 2"
    t_datas3 = calculate_for_model(m, ThinDisc(0.0, Inf), LampPostModel(h=5.0))
    @info "Thin disc 3"
    t_datas4 = calculate_for_model(m, ThinDisc(0.0, Inf), LampPostModel(h=10.0))
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
    x_lims1 = (0.85, 1.05)
    x_lims2 = (0.65, 1.2)
    x_lims3 = (0.45, 1.4)
    x_lims4 = (0.45, 1.6)

    axes = map(enumerate((m_datas2, m_datas3, m_datas4))) do (row, ds)
        map(enumerate(ds)) do (col, d)
            ax = if col == 1
                Axis(ga[row, col],
                    ylabel ="τ",
                    xlabel ="E/E₀",
                    yticks = LinearTicks(3),
                )
            else
                Axis(ga[row, col],
                    xlabel ="E/E₀",
                    yticks = LinearTicks(3),
                )
            end
            thin = (t_datas2, t_datas3, t_datas4)[row][col]

            for (lim) in ((lims1, lims2, lims3, lims4))
                thin_le = lag_energy(thin.all_pulses, lim...)
                lines!(ax, thin.E, thin_le, color = :grey, alpha=0.3)

                le = lag_energy(d.all_pulses, lim...)
                lines!(ax, d.E, le)
            end
            if row != 3 
                hidexdecorations!(ax, grid=false)
            end
            xlims!(ax, (x_lims1, x_lims2, x_lims3, x_lims4)[col]...)
            ax
        end
    end
    axs = reshape(reduce(vcat, axes), (length(axes[1]), length(axes)))
    for i in 1:length(axes[1])
        linkxaxes!(axs[i, 1], axs[i, 2])
        # linkxaxes!(axs[3, i], axs[4, i])
    end
    # for i in 1:length(axes)
    #     root = axs[1, i]
    #     for a in axs[2:end, i]
    #         linkyaxes!(root, a)
    #     end
    # end
    
    rowgap!(ga, 10)
    colgap!(ga, 10)

    for (i, a) in enumerate([10, 30, 60, 80])
        Label(
            ga[1, i, Top()],
            text = "θ = $(a)°",
            padding = (0, 0, 10, 0),
            fontsize = 13,
            font = :bold,
        )
    end
    for (i, h) in enumerate([2, 5, 10])
        Label(
            ga[i, 4, Right()],
            text = "h = $h",
            padding = (10, 0, 00, 0),
            fontsize = 13,
            font = :bold,
        )
    end

    @savefigure(fig)
    fig
end
