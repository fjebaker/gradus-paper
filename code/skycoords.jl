using Makie, CairoMakie, LaTeXStrings, ColorSchemes, StaticArrays

include("common.jl")
include("makie-drawing.jl")
projector(θ, ϕ) = LinearMap(RotZ(ϕ)) ∘ LinearMap(RotX(θ))

abstract type AbstractDrawable end

struct SkySphere <: AbstractDrawable
    r::Float64
end
draw(ax, s::SkySphere; kwargs...) = bounding_sphere!(ax, R = s.r; kwargs...)

struct SkyRay <: AbstractDrawable
    θ::Float64
    ϕ::Float64
    l::Float64
end

endpoint(r::SkyRay) = projector(r.θ, r.ϕ)([0, r.l, 0])

function draw(ax, r::SkyRay; kwargs...)
    p = endpoint(r)
    ps = hcat(zeros(eltype(p), 3), p)
    lines!(ax, ps[1, :], ps[2, :], ps[3, :]; kwargs...)
end

struct PointOnSphere <: AbstractDrawable
    s::SkySphere
    θ::Float64
    ϕ::Float64
end
draw(ax, p::PointOnSphere; kwargs...) =
    scatter!(ax, projector(p.θ, p.ϕ)([0, p.s.r, 0]); kwargs...)

abstract type AbstractOrientation end
struct Horizontal <: AbstractOrientation end
struct Vertical <: AbstractOrientation end

struct SphereArc{O<:AbstractOrientation} <: AbstractDrawable
    s::SkySphere
    orientation::O
    θ::Float64
    ϕ::Float64
end

get_projection(arc::SphereArc{<:Horizontal}) =
    LinearMap(RotZ(arc.ϕ)) ∘ LinearMap(RotX(arc.θ))
get_projection(arc::SphereArc{<:Vertical}) = LinearMap(RotZ(arc.ϕ)) ∘ LinearMap(RotY(π / 2))

function draw(ax, arc::SphereArc; kwargs...)
    t = get_projection(arc)
    R = arc.s.r
    phi_circ = 0.0:0.001:2π
    x = @. R * cos(phi_circ)
    y = @. R * sin(phi_circ)
    z = zeros(length(y))

    points = reduce(hcat, [x, y, z])'
    translated = reduce(hcat, map(t, eachcol(points)))

    plot_line_occluded!(
        ax,
        translated[1, :],
        translated[2, :],
        translated[3, :],
        R;
        linewidth = 1.9,
        kwargs...,
    )
end

struct ImagePlane <: AbstractDrawable
    s::SkySphere
    θ::Float64
    ϕ::Float64
    x_dim::Float64
    z_dim::Float64
end

function midpoint(plane::ImagePlane)
    T = projector(plane.θ, plane.ϕ)
    T(SVector(0, plane.s.r, 0))
end
function draw(ax, plane::ImagePlane; kwargs...)
    T = projector(plane.θ, plane.ϕ)
    ps = map(
        T,
        (
            [plane.x_dim, plane.s.r, plane.z_dim],
            [-plane.x_dim, plane.s.r, plane.z_dim],
            [-plane.x_dim, plane.s.r, -plane.z_dim],
            [plane.x_dim, plane.s.r, -plane.z_dim],
        ),
    )
    ps = (ps..., ps[1])
    p = reduce(hcat, ps)
    lines!(ax, p[1, :], p[2, :], p[3, :]; kwargs...)
end

struct CoordinateLines{T} <: AbstractDrawable
    ob::T
end
coordlines(im::ImagePlane) = CoordinateLines(im)
function draw(ax, cl::CoordinateLines{<:ImagePlane}; kwargs...)
    T = projector(cl.ob.θ, cl.ob.ϕ)
    X = reduce(hcat, map(T, ([-cl.ob.x_dim, cl.ob.s.r, 0], [cl.ob.x_dim, cl.ob.s.r, 0])))
    Z = reduce(hcat, map(T, ([0, cl.ob.s.r, -cl.ob.z_dim], [0, cl.ob.s.r, cl.ob.z_dim])))

    lines!(ax, X[1, :], X[2, :], X[3, :]; kwargs...)
    lines!(ax, Z[1, :], Z[2, :], Z[3, :]; kwargs...)
end

struct PointOnImage <: AbstractDrawable
    p::ImagePlane
    θ::Float64
    ϕ::Float64
end

function tangent_project(n, p, o, d)
    t = ((p .- o) ⋅ n) / (d ⋅ n)
    t .* d
end

function point_on_plane(p::PointOnImage)
    d = spher_to_cart(1, -p.θ + π / 2, p.ϕ + π / 2)
    n = spher_to_cart(1, -p.p.θ + π / 2, p.p.ϕ + π / 2)
    p = spher2cart(p.p.s.r, -p.p.θ + π / 2, p.p.ϕ + π / 2)
    tangent_project(n, p, zeros(Float64, 3), d)
end

function draw(ax, p::PointOnImage; kwargs...)
    x = point_on_plane(p)
    scatter!(ax, x; kwargs...)
end

intersection(s::SkySphere, r::SkyRay) = PointOnSphere(s, r.θ, r.ϕ)
intersection(s::ImagePlane, r::SkyRay) = PointOnImage(s, r.θ, r.ϕ)

struct Orthogonal{T} <: AbstractDrawable
    l1::T
    l2::T
    s::Float64
end
function draw(ax, orth::Orthogonal{<:SkyRay}; kwargs...)
    # get unit vectors
    p1 = projector(orth.l1.θ, orth.l1.ϕ)([0, orth.s, 0])
    p2 = projector(orth.l2.θ, orth.l2.ϕ)([0, orth.s, 0])
    x = p1 .+ p2
    ps = hcat(p1, x, p2)
    lines!(ax, ps[1, :], ps[2, :], ps[3, :]; kwargs...)
end

to_unit(v) = v ./ √(v ⋅ v)

struct FiniteArc{V} <: AbstractDrawable
    v::V
    w::V
end

FiniteArc(v1, v2, l) = FiniteArc(to_unit(v1) .* l, to_unit(v2) .* l)

to_vector(a::SkyRay) = projector(a.θ, a.ϕ)([0, a.l, 0])

function FiniteArc(θ1, ϕ1, θ2, ϕ2, l)
    v = projector(θ1, ϕ1)([0, l, 0])
    w = projector(θ2, ϕ2)([0, l, 0])
    FiniteArc(v, w)
end

function draw(ax, f::FiniteArc; N = 50, kwargs...)
    v = f.v
    w = f.w
    v2 = v ⋅ v
    w2 = w ⋅ w
    angdiff = acos((v ⋅ w) / √(v2 * w2))
    n = to_unit(v × w) .* (angdiff / N)

    T = LinearMap(RotationVec(n[1], n[2], n[3]))
    points = [v]
    for _ = 1:N
        push!(points, T(points[end]))
    end
    ps = reduce(hcat, points)

    lines!(ax, ps[1, :], ps[2, :], ps[3, :]; kwargs...)
end

function draw(ax, x::SVector{3}; kwargs...)
    scatter!(ax, x; kwargs...)
end

struct PointLine{P} <: AbstractDrawable
    points::P
end
function PointLine(p1, p2, ps...)
    PointLine((p1, p2, ps...))
end

function draw(ax, p::PointLine; kwargs...)
    ps = reduce(hcat, p.points)
    lines!(ax, ps[1, :], ps[2, :], ps[3, :]; kwargs...)
end

function axes_points(k::PointOnImage)
    T = projector(k.p.θ, k.p.ϕ)
    a1 = to_unit(T(SVector(k.p.x_dim, 0, 0)))
    a2 = to_unit(T(SVector(0, 0, k.p.z_dim)))
    O = midpoint(k.p)
    x = point_on_plane(k)
    A = (a1 ⋅ x) .* a1 .+ O
    B = (a2 ⋅ x) .* a2 .+ O
    (A, B, SVector(x))
end

function point_between(p1::SVector{3}, p2::SVector{3}; t = 0.5)
    @. t * (p2 - p1) + p1
end

angle_between(a::SkyRay, b::SkyRay, l) = FiniteArc(a.θ, a.ϕ, b.θ, b.ϕ, l)

begin
    # palette = Iterators.Stateful(Iterators.Cycle(ColorSchemes.tab10.colors))
    palette = Iterators.Stateful(Iterators.Cycle(Makie.wong_colors()))

    dim = 15
    fig = Figure(resolution = (650, 450))
    ax = Makie.Axis3(
        fig[1, 1],
        # aspect = (1.0, 1.0, 0.6),
        aspect = :data,
        limits = (-dim, dim, -dim, dim, -0.30dim, 0.8dim),
        elevation = π / 20,
        # elevation = π / 6,
        azimuth = deg2rad(90 + 20),
        xspinewidth = 0,
        yspinewidth = 0,
        zspinewidth = 0,
        # salignmode=Outside(),
        # protrusions = 0,
        viewmode = :stretch,
        # backgroundcolor = RGBA(0.0,0.0,0.0,0.0),
        xgridvisible = false,
        ygridvisible = false,
        zgridvisible = false,
        xlabelvisible = false,
        ylabelvisible = false,
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0, 0),
        zautolimitmargin = (0, 0),
        zticklabelpad = 0,
    )
    Makie.hidedecorations!(ax)

    # lines!(ax, [0, 0], [0, 10], [0, 0], color = :red)

    phi = deg2rad(42)
    theta = deg2rad(67)

    color_sphere = popfirst!(palette)
    color_plane = popfirst!(palette)
    color_tetrad = popfirst!(palette)
    # _ = popfirst!(palette)
    color_ray = popfirst!(palette)
    # _ = popfirst!(palette)
    color_helpers = popfirst!(palette)
    _ = popfirst!(palette)
    _ = popfirst!(palette)
    color_impact = popfirst!(palette)
    # color_upsilon = :red
    color_upsilon = popfirst!(palette)
    _ = popfirst!(palette)
    color_Psi = popfirst!(palette)

    # color_sphere = :black
    # color_plane = :black
    # color_tetrad = :black
    # color_ray = :black
    # color_helpers = :black
    # color_impact = :black
    # color_upsilon = :black
    # color_Psi = :black

    sphere = SkySphere(10)
    draw(ax, sphere, color = color_sphere, linewidth = 3.0)

    p1 = PointOnSphere(sphere, theta, phi)

    plane = ImagePlane(sphere, theta, phi, 14.8, 9.8)
    draw(ax, plane, color = color_plane, linewidth = 2.0)
    plane_coords = coordlines(plane)
    draw(ax, plane_coords, color = color_plane, linestyle = :dash)

    arc1 = SphereArc(sphere, Horizontal(), theta, phi)
    draw(ax, arc1; color = color_sphere, linewidth = 1.0)
    arc2 = SphereArc(sphere, Vertical(), theta, phi)
    draw(ax, arc2; linewidth = 1.0, color = color_sphere)

    baseline = SphereArc(sphere, Horizontal(), theta - π / 2, phi)
    draw(ax, baseline; linewidth = 1.0, linestyle = :dash, color = color_sphere)

    prad = SkyRay(theta, phi, 13)
    draw(ax, prad, color = color_tetrad)

    ray2 = SkyRay(theta - deg2rad(40), phi + deg2rad(47), 21)
    draw(ax, ray2, color = color_ray)

    draw(
        ax,
        intersection(sphere, ray2),
        color = RGBAf(0.0, 0.0, 0.0, 0.0),
        strokewidth = 1.0,
        strokecolor = color_ray,
    )
    pinter = intersection(plane, ray2)

    k1, k2, P = axes_points(pinter)
    draw(ax, PointLine(k1, P, k2); color = color_impact)
    O = midpoint(plane)

    draw(ax, PointLine(SVector(0.0, 0.0, 0.0), P .- O), color = color_helpers)
    draw(
        ax,
        PointLine(P .- O, P),
        linewidth = 1.0,
        linestyle = :dash,
        color = color_helpers,
    )

    p_plane = point_on_plane(pinter)

    ptheta = SkyRay(theta - deg2rad(90), phi, prad.l)
    draw(ax, ptheta, color = color_tetrad)
    pphi = SkyRay(0, phi + deg2rad(90), prad.l)
    draw(ax, pphi, color = color_tetrad)

    draw(ax, Orthogonal(ptheta, pphi, 0.8); color = color_tetrad)
    draw(ax, Orthogonal(ptheta, prad, 0.8); color = color_tetrad)
    draw(ax, Orthogonal(pphi, prad, 0.8); color = color_tetrad)

    upsilon_arc = angle_between(prad, ray2, 3.0)
    draw(ax, upsilon_arc; color = color_upsilon)

    Psi = FiniteArc(to_vector(pphi), P .- O, 8.0)
    draw(ax, Psi; color = color_Psi)

    draw(ax, pinter, color = color_ray, marker = :x)

    text!(ax, point_between(k2, P), text = L"\alpha", fontsize = 20, color = color_impact)
    text!(ax, point_between(k1, P), text = L"\beta", fontsize = 20, color = color_impact)
    text!(ax, endpoint(prad), text = L"e_{(r)}", fontsize = 20, color = :black)
    text!(ax, endpoint(ptheta), text = L"e_{(\theta)}", fontsize = 20, color = :black)
    text!(ax, endpoint(pphi), text = L"e_{(\phi)}", fontsize = 20, color = :black)

    text!(
        ax,
        point_between(upsilon_arc.w, upsilon_arc.v),
        text = L"\Upsilon",
        fontsize = 20,
        color = color_upsilon,
    )
    text!(
        ax,
        point_between(Psi.v, Psi.w, t = 0.8),
        text = L"\Psi",
        fontsize = 20,
        color = color_Psi,
        align = (:right, :bottom),
    )
    text!(ax, endpoint(ray2), text = L"v_{(\mu)}", fontsize = 20, color = color_ray)


    resize_to_layout!(fig)
    @savefigure(fig)
    fig
end
