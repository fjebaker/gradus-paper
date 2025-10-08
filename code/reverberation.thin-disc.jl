using Gradus, Makie, CairoMakie, LaTeXStrings, Printf

include("common.jl")

function _format_model(model)
    hh = Printf.@sprintf "%.0f" model.h
    L"h = %$hh r_\text{g}"
end

function classical_time(theta, r, h)
    v = r^2 - 2r * h * cos(theta) + h^2
    sqrt(v)
end

function calculate_2d_transfer_function(m, x, model, itb, prof, radii)
    bins = collect(range(0.0, 1.5, 300))
    tbins = collect(range(0, 2000.0, 3000))

    t0 = continuum_time(m, x, model)
    @show t0

    flux = @time Gradus.integrate_lagtransfer(
        prof,
        itb,
        bins,
        tbins;
        t0 = t0,
        n_radii = 10000,
        rmin = minimum(radii),
        rmax = maximum(radii),
        h = 1e-8,
        g_grid_upscale = 10,
    )

    flux[flux.==0] .= NaN
    bins, tbins, flux
end

function calculate_lag_transfer(m, d, model, radii, itb)
    prof = @time emissivity_profile(m, d, model; n_samples = 100_000)
    E, t, f = @time calculate_2d_transfer_function(m, x, model, itb, prof, radii)
    ψ = Gradus.sum_impulse_response(f)
    freq, τ = @time lag_frequency(t, f)
    freq, τ, ψ, t
end

m = KerrMetric(1.0, 0.998)
x = SVector(0.0, 10_000.0, deg2rad(45), 0.0)
radii = Gradus.Grids._inverse_grid(Gradus.isco(m), 1000.0, 400)

# models
model1 = LampPostModel(h = 2.0)
model2 = LampPostModel(h = 5.0)
model3 = LampPostModel(h = 10.0)
model4 = LampPostModel(h = 20.0)

# thin disc
d = ThinDisc(0.0, Inf)

itb = Gradus.interpolated_transfer_branches(m, x, d, radii; verbose = true)

freq1, τ1, impulse1, time1 = calculate_lag_transfer(m, d, model1, radii, itb)
freq2, τ2, impulse2, time2 = calculate_lag_transfer(m, d, model2, radii, itb)
freq3, τ3, impulse3, time3 = calculate_lag_transfer(m, d, model3, radii, itb)
freq4, τ4, impulse4, time4 = calculate_lag_transfer(m, d, model4, radii, itb)

# thick disc
thick_d = ShakuraSunyaev(m)

thick_itb = Gradus.interpolated_transfer_branches(m, x, d, radii; verbose = true, β₀ = 2.0)

thick_freq1, thick_τ1, thick_impulse1, thick_time1 =
    calculate_lag_transfer(m, thick_d, model1, radii, thick_itb)
thick_freq2, thick_τ2, thick_impulse2, thick_time2 =
    calculate_lag_transfer(m, thick_d, model2, radii, thick_itb)
thick_freq3, thick_τ3, thick_impulse3, thick_time3 =
    calculate_lag_transfer(m, thick_d, model3, radii, thick_itb)
thick_freq4, thick_τ4, thick_impulse4, thick_time4 =
    calculate_lag_transfer(m, thick_d, model4, radii, thick_itb)

begin
    @show sum(filter(!isnan, impulse1))
    @show sum(filter(!isnan, impulse2))
    @show sum(filter(!isnan, impulse3))
    @show sum(filter(!isnan, impulse4))
end

begin
    palette = _default_palette()

    fig = Figure(resolution = (500, 650))
    ga = fig[1, 1] = GridLayout()
    ax1 =
        Axis(ga[2, 1], yscale = log10, xlabel = L"Time $t$", ylabel = L"Impulse Response$$")
    ax2 = Axis(
        ga[3, 1],
        xscale = log10,
        xminorgridvisible = true,
        xminorticks = IntervalsBetween(10),
        ylabel = L"Lag $\tau$",
        xlabel = L"Frequency $f$",
    )

    color = popfirst!(palette)
    lines!(ax1, time1, abs.(impulse1), color = color)
    lines!(ax1, thick_time1, abs.(thick_impulse1), color = color, linestyle = :dot)
    l1 = lines!(ax2, freq1, τ1, label = _format_model(model1), color = color)
    lines!(
        ax2,
        thick_freq1,
        thick_τ1,
        label = _format_model(model1),
        linestyle = :dot,
        color = color,
    )

    color = popfirst!(palette)
    lines!(ax1, time2, abs.(impulse2), color = color)
    lines!(ax1, thick_time2, abs.(thick_impulse2), color = color, linestyle = :dot)
    l2 = lines!(ax2, freq2, τ2, label = _format_model(model2), color = color)
    lines!(
        ax2,
        thick_freq2,
        thick_τ2,
        label = _format_model(model2),
        linestyle = :dot,
        color = color,
    )

    color = popfirst!(palette)
    lines!(ax1, time3, abs.(impulse3), color = color)
    lines!(ax1, thick_time3, abs.(thick_impulse3), color = color, linestyle = :dot)
    l3 = lines!(ax2, freq3, τ3, label = _format_model(model3), color = color)
    lines!(
        ax2,
        thick_freq3,
        thick_τ3,
        label = _format_model(model3),
        linestyle = :dot,
        color = color,
    )

    color = popfirst!(palette)
    lines!(ax1, time4, abs.(impulse4), color = color)
    lines!(ax1, thick_time4, abs.(thick_impulse4), color = color, linestyle = :dot)
    l4 = lines!(ax2, freq4, τ4, label = _format_model(model4), color = color)
    lines!(
        ax2,
        thick_freq4,
        thick_τ4,
        label = _format_model(model4),
        linestyle = :dot,
        color = color,
    )

    lines!(ax2, [1e-5, 1], [0, 0], color = :black)

    xx = collect(range(1e-5, 1, 1000))
    yy = @. 1 / (2 * π * x)

    Legend(
        ga[1, 1],
        [l1, l2, l3, l4],
        map(_format_model, [model1, model2, model3, model4]),
        orientation = :horizontal,
        height = 10,
        framevisible = false,
        padding = (0, 0, 0, 0),
    )

    xlims!(ax2, 5e-5, 0.3)
    ylims!(ax2, -10, 50)

    ylims!(ax1, 2e-5, 0.2)
    xlims!(ax1, 0, 120)

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

    @savefigure(fig)
    fig
end
