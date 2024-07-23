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
    δr = 2.0
    if r ≤ R
        one(T)
    elseif R ≤ r ≤ R + δr
        t = (r - R) / δr
        # use an arbitrarily steep smooth interpolation
        # this one isn't perfect, but does a good job
        k = atan(1e5t) * 2 / π
        (1 - k)
    else
        zero(T)
    end
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
    k = transition_function(r, m.R)
    if k >= 1
        Gradus.metric_components(m.m1, rθ)
    elseif k <= 0
        Gradus.metric_components(m.m2, rθ)
    else
        c1 = Gradus.metric_components(m.m1, rθ)
        c2 = Gradus.metric_components(m.m1, rθ)
        @. c1 * k + (1 - k) * c2
    end
end
Gradus.inner_radius(m::SwitchingMetric) = Gradus.inner_radius(m.m1)
Gradus.isco(m::SwitchingMetric) = Gradus.isco(m.m1)

begin
    heights = Gradus.Grids._geometric_grid(2.0, 100.0, 200) |> collect
    push!(heights, 100.0)
end

x = SVector(0.0, 10_000.0, deg2rad(45), 0.0)

s_m1 = SwitchingMetric(KerrMetric(1.0, 0.998), Gradus.SphericalMetric(), 100.0)
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
    l1 = lines!(ax, heights, abs.(s_times1 .- s_times0), color = c1)
    l2 = lines!(ax, heights, abs.(s_times2 .- s_times0), color = c2)
    l3 = lines!(ax, heights, abs.(s_times3 .- s_times0), color = c3)

    Legend(
        ga[1, 1],
        [l3, l2, l1],
        [L"R=20", L"R=50", L"R=100"],
        orientation = :horizontal,
        padding = (0, 0, 0, 0),
        framevisible = false,
    )

    rowgap!(ga, 10)

    @savefigure(fig)
    fig
end
