using Gradus, Makie, CairoMakie, LaTeXStrings, Printf

include("common.jl")

function _format_angle(x::SVector{4})
    fmt = Printf.@sprintf "%.0f" rad2deg(x[3])
    L"%$fmt^\circ"
end

function calculate_times(m, heights, x)
    @time map(heights) do h
        try
            continuum_time(m, x, LampPostModel(h = h, θ = 0.001))
        catch
            NaN
        end
    end
end

function calculate_reflection_time(m, heights, x)
    d = ThinDisc(0.0, Inf)
    radii = collect(range(Gradus.isco(m), 50.0, 50))
    itb = Gradus.interpolated_transfer_branches(m, x, d, radii, verbose = true)

    bins = collect(range(0.0, 1.5, 50))
    tbins = collect(range(0, 100.0, 500))

    map(heights) do h
        prof = emissivity_profile(m, d, LampPostModel(h = h))
        flux = @time Gradus.integrate_lagtransfer(
            prof,
            itb,
            radii,
            bins,
            tbins;
            t0 = x[2],
            Nr = 500,
            h = 1e-8,
        )
        t = Gradus.sum_impulse_response(flux)
        i = findfirst(>(0), t)
        tbins[i] + x[2]
    end
end

function transition_function(r::T, R) where {T}
    Gradus._smooth_interpolate(r, R; δx = 10.0)
end

struct SwitchingMetric{M1,M2,T} <: AbstractStaticAxisSymmetric{T}
    m1::M1
    m2::M2
    R::T
    function SwitchingMetric(m1::AbstractMetric{T}, m2::AbstractMetric{T}, R::T) where {T}
        new{typeof(m1),typeof(m2),T}(m1, m2, R)
    end
end

function Gradus.metric_components(m::SwitchingMetric, rθ)
    r = rθ[1]
    c1 = Gradus.metric_components(m.m1, rθ)
    c2 = Gradus.metric_components(m.m2, rθ)
    k = clamp(transition_function(r, m.R), 0, 1)
    @. c1 * k + (1 - k) * c2
end
Gradus.inner_radius(m::SwitchingMetric) = Gradus.inner_radius(m.m1)
Gradus.isco(m::SwitchingMetric) = Gradus.isco(m.m1)

begin
    heights = range(2.0, 110.0, 200) |> collect
end

x = SVector(0.0, 10_000.0, deg2rad(45), 0.0)

s_m1 = SwitchingMetric(KerrMetric(1.0, 0.998), KerrMetric(1e-5, 0.0), 100.0)
s_m2 = SwitchingMetric(s_m1.m1, s_m1.m2, 50.0)
s_m3 = SwitchingMetric(s_m1.m1, s_m1.m2, 20.0)

s_times0 = calculate_times(s_m1.m1, heights, x)
s_times1 = calculate_times(s_m1, heights, x)
s_times2 = calculate_times(s_m2, heights, x)
s_times3 = calculate_times(s_m3, heights, x)

begin
    fig = Figure(resolution = (450, 350))
    ga = fig[1, 1] = GridLayout()
    ax = Axis(
        ga[2, 1],
        xlabel = L"h",
        ylabel = L"\delta t",
        xscale = log10,
        xticks = [1, 2, 5, 50, 10, 100],
        xminorticks = [3, 4, 6, 7, 8, 9, 20, 30, 40, 60, 70, 80, 90],
        xminorgridvisible = true,
    )

    palette = _default_palette()
    c1 = popfirst!(palette)
    c2 = popfirst!(palette)
    c3 = popfirst!(palette)
    l3 = lines!(ax, heights, (s_times3 .- s_times0), color = c3)
    l2 = lines!(ax, heights, (s_times2 .- s_times0), color = c2)
    l1 = lines!(ax, heights, (s_times1 .- s_times0), color = c1)

    Legend(
        ga[1, 1],
        [l3, l2, l1],
        [L"R=20", L"R=50", L"R=100"],
        orientation = :horizontal,
        padding = (0, 0, 0, 0),
        framevisible = false,
    )

    rowgap!(ga, 10)

    # @savefigure(fig)
    fig
end


function _format_model(model)
    hh = Printf.@sprintf "%.0f" model.h
    L"h = %$hh r_\text{g}"
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
        n_radii = 8000,
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
radii = Gradus.Grids._inverse_grid(Gradus.isco(m), 1000.0, 200)

# models
model2 = LampPostModel(h = 10.0)

# thin disc
d = ThinDisc(0.0, Inf)

itb = Gradus.interpolated_transfer_branches(m, x, d, radii; verbose = true)

freq2, τ2, impulse2, time2 = calculate_lag_transfer(m, d, model2, radii, itb)

# alternating metric
alt_m = s_m1

Gradus._redshift_guard(m::SwitchingMetric, args...) = Gradus._redshift_guard(m.m1, args...)
Gradus.ConstPointFunctions.redshift(m::SwitchingMetric, args...) = Gradus.ConstPointFunctions.redshift(m.m1, args...)

# gridding for the photon plane
plane = PolarPlane(GeometricGrid(); Nr = 1200, Nθ = 1200, r_max = 700.0)

# integrate source to disc and observer to disc
tf = @time lagtransfer(
    alt_m,
    x,
    ThinDisc(Gradus.isco(m), Inf),
    model2
    ;
    plane = plane,
    callback = domain_upper_hemisphere(),
    n_samples = 1,
    verbose = true,
)
prof = @time emissivity_profile(alt_m, ThinDisc(Gradus.isco(m), Inf), model2; n_samples = 1_000)
alt_t0 = continuum_time(alt_m, x, model2)

begin
    alt_time, alt_E, alt_f = binflux(tf, prof; N_E = 1500, N_t = 3000, t0 = alt_t0, time_lims = (alt_t0, alt_t0 + 5000.0))
    replace!(alt_f, NaN => 0.0)
    Gradus._normalize!(alt_f, alt_E)
    alt_impulse = Gradus.sum_impulse_response(alt_f)
    alt_freq, alt_τ = @time lag_frequency(alt_time, alt_f)
end

begin
    palette = _default_palette()

    fig = Figure(resolution = (500, 450))
    ga = fig[1, 1] = GridLayout()
    ax2 = Axis(
        ga[1, 1],
        xscale = log10,
        xminorgridvisible = true,
        xminorticks = IntervalsBetween(10),
        ylabel = L"Lag $\tau$",
        xlabel = L"Frequency $f$",
    )

    color = popfirst!(palette)
    l2 = lines!(ax2, freq2, τ2, label = _format_model(model2), color = color)
    lines!(
        ax2,
        alt_freq,
        alt_τ,
        label = _format_model(model2),
        linestyle = :dot,
        color = color,
    )

    lines!(ax2, [1e-5, 1], [0, 0], color = :black)

    xx = collect(range(1e-5, 1, 1000))
    yy = @. 1 / (2 * π * x)

    xlims!(ax2, 5e-5, 0.3)
    ylims!(ax2, -10, 50)

    # @savefigure(fig)
    fig
end