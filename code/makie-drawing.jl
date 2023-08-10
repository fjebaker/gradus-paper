using LinearAlgebra
import Makie
import CairoMakie
using CoordinateTransformations, Rotations

function clip_plot_lines!(ax, x, y, z; dim = 10, kwargs...)
    mask = @. (x > dim) | (x < -dim) | (y > dim) | (y < -dim) | (z > 1.05dim)
    mask = .!mask
    Makie.lines!(ax, x[mask], y[mask], z[mask]; kwargs...)
end
function clip_plot_scatter!(ax, x, y, z; R = 1, kwargs...)
    if first(is_visible(ax, [x; y; z], R))
        clip_plot_scatter!(ax, [x], [y], [z]; kwargs...)
    end
end
function clip_plot_scatter!(
    ax,
    x::AbstractArray,
    y::AbstractArray,
    z::AbstractArray;
    dim = 10,
    kwargs...,
)
    mask = @. (x > dim) | (x < -dim) | (y > dim) | (y < -dim) | (z > 1.1dim)
    mask = .!mask
    Makie.scatter!(ax, x[mask], y[mask], z[mask]; kwargs...)
end

spher_to_cart(r, θ, ϕ) = (r * sin(θ) * cos(ϕ), r * sin(θ) * sin(ϕ), r * cos(θ))

function plot_obscured_line!(ax, A, B; R = 1, kwargs...)
    p = A
    dvec = B - A
    for i = 1:1000
        if first(is_visible(ax, [p[1]; p[2]; p[3]], R)) && (LinearAlgebra.norm(p) > R)
            break
        end
        p += dvec * 0.001
    end

    x, y, z = [p[1], B[1]], [p[2], B[2]], [p[3], B[3]]
    clip_plot_lines!(ax, x, y, z; linewidth = 0.6, kwargs...)
end

function plot_line_occluded!(ax, x, y, z, R; kwargs...)
    points_ = reduce(hcat, [x, y, z])'
    mask = is_visible(ax, points_, R)

    s = 1
    e = 1
    for i in mask
        if !i
            if s != e
                clip_plot_lines!(
                    ax,
                    x[s:e-1],
                    y[s:e-1],
                    z[s:e-1];
                    linewidth = 0.8,
                    kwargs...,
                )
            end
            s = e
        end
        e += 1
    end
    if s != e
        clip_plot_lines!(ax, x[s:e-1], y[s:e-1], z[s:e-1]; linewidth = 0.8, kwargs...)
    end
end

function circle_points(θ; N = 100, r = 1.0)
    ϕ = range(0.0, 1, N) * 2π # 0.0:0.01:2π
    x = @. r * sin(θ) * cos(ϕ)
    y = @. r * sin(θ) * sin(ϕ)
    z = r * cos(θ) .* ones(length(ϕ))
    x, y, z
end

function bounding_sphere!(ax; R = 1.005, kwargs...)
    phi_circ = 0.0:0.001:2π
    x = @. R * cos(phi_circ)
    y = @. R * sin(phi_circ)
    z = zeros(length(y))

    t = LinearMap(RotZ(ax.azimuth.val - π)) ∘ LinearMap(RotY(ax.elevation.val - π / 2))
    points = reduce(hcat, [x, y, z])'
    translated = reduce(hcat, map(t, eachcol(points)))


    clip_plot_lines!(
        ax,
        translated[1, :],
        translated[2, :],
        translated[3, :];
        linewidth = 1.9,
        kwargs...,
    )
    #translated
end

function is_visible(ax, points, R)
    # transform viewing angle to normal vector in data coordinates
    a = ax.azimuth.val - π
    e = ax.elevation.val - π / 2
    n = [spher_to_cart(R, e, a)...]
    n = n ./ LinearAlgebra.norm(n)
    # clip_plot_scatter!(ax, [n[1]], [n[2]], [n[3]])
    map(eachcol(points)) do p
        k = (p ⋅ n)
        if k > 0 # infront
            true
        else
            Pp = p - (k .* n)
            sqrt(abs(Pp ⋅ Pp)) > R
        end
    end
end

function plot_sol(ax, x; R = 1.0, kwargs...)
    cart_points = []
    for t in range(x.t[1], x.t[end], 2_000)
        p = x(t)[1:4]
        if p[2] > (1.08 * R)
            push!(cart_points, [spher2cart(p...)...])
        end
    end
    res = reduce(hcat, cart_points)'
    # Plots.plot3d!(res[:, 1], res[:, 2], res[:, 3] ; color = COLOR)
    # plot_line_occluded!(ax, res[:, 1], res[:, 2], res[:, 3], R)
    mask = is_visible(ax, res', R)

    # plot disconnected continuous regions
    s = 1
    e = findnext(==(0), mask, s)
    if isnothing(e)
        e = lastindex(mask)
    end
    @views while e < lastindex(mask)
        clip_plot_lines!(ax, res[s:e, 1], res[s:e, 2], res[s:e, 3], ; kwargs...)
        s = findnext(==(true), mask, e)
        if isnothing(s)
            s = lastindex(mask)
        end
        e = findnext(==(false), mask, s)
        if isnothing(e)
            e = lastindex(mask)
        end
    end
    clip_plot_lines!(ax, res[s:e, 1], res[s:e, 2], res[s:e, 3], ; kwargs...)
    if unpack_solution(x).status == StatusCodes.IntersectedWithGeometry
        clip_plot_scatter!(
            ax,
            res[end, 1],
            res[end, 2],
            res[end, 3],
            ;
            R = R,
            markersize = 6,
            kwargs...,
        )
    end
end

function plotring(ax, r, z0; R, kwargs...)
    x = Float64[]
    y = Float64[]
    z = Float64[]
    for ϕ in range(0.0, 2π, 500)
        X, Y, Z = spher2cart(r, π / 2, ϕ)
        Z += z0
        push!(x, X)
        push!(y, Y)
        push!(z, Z)
    end
    plot_line_occluded!(ax, x, y, z, R; linewidth = 0.8, kwargs...)
end

spher2cart(r, θ, ϕ) = r * sin(θ) * cos(ϕ), r * sin(θ) * sin(ϕ), r * cos(θ)
cart2sphere(x, y, z) = (√(x^2 + y^2 + z^2), atan(y, x), atan(√(x^2 + y^2), z))
spher2cart(_, r, θ, ϕ) = spher2cart(r, θ, ϕ)
