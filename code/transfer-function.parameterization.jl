using Makie, CairoMakie, LaTeXStrings, Printf
using Gradus

include("common.jl")

SMALL_FONT_SIZE = 18
LARGE_FONT_SIZE = 24

struct RedshiftBranch{K}
    α_interp::K
    β_interp::K
end

function _make_branch(g✶::AbstractArray, ps::AbstractArray)
    I = sortperm(g✶)
    α = first.(ps)[I]
    β = last.(ps)[I]
    t1 = Gradus.DataInterpolations.LinearInterpolation(α, g✶[I])
    t2 = Gradus.DataInterpolations.LinearInterpolation(β, g✶[I])
    RedshiftBranch(t1, t2)
end

struct RedshiftData{T,K}
    visible::BitVector
    α::Vector{T}
    β::Vector{T}
    upper::RedshiftBranch{K}
    lower::RedshiftBranch{K}
    gmin::T
    gmax::T
    r::T
end

function RedshiftData(α, β, g, r, I)
    gmin, imin = findmin(g)
    gmax, imax = findmax(g)

    g✶ = Gradus.g_to_g✶.(g, gmin, gmax)
    points = [SVector(a, b) for (a, b) in zip(α, β)]
    gup, pup, gdown, pdown = splitbranches(g✶, points, imin, imax)

    A = _make_branch(gup, pup)
    B = _make_branch(gdown, pdown)

    RedshiftData(I, α, β, A, B, gmin, gmax, r)
end

function interpolate(rd::RedshiftData, g✶, upper)
    rb = (upper ? rd.upper : rd.lower)
    α, β = rb.α_interp(g✶), rb.β_interp(g✶)
    # work out if it's visible by comparing it to closest known α, β
    _, I = findmin(i -> (rd.α[i] - α)^2 + (rd.β[i] - β)^2, eachindex(rd.α))
    rd.visible[I] ? SVector(α, β) : SVector(NaN, NaN)
end

function splitbranches(g✶::Vector{T}, f::Vector{V}, imin, imax) where {T,V}
    i1, i2 = imax > imin ? (imin, imax) : (imax, imin)
    if (i1 == i2)
        error("Resolved same min/max")
    end

    N1 = i2 - i1 + 1
    N2 = length(f) - N1 + 2
    branch1_f = zeros(V, N1)
    branch2_f = zeros(V, N2)
    branch1_g✶ = zeros(Float64, N1)
    branch2_g✶ = zeros(Float64, N2)

    for (i, j) in enumerate(i1:i2)
        branch1_f[i] = f[j]
        branch1_g✶[i] = g✶[j]
    end
    for (i, j) in enumerate(Iterators.flatten((1:i1, i2:length(f))))
        branch2_f[i] = f[j]
        branch2_g✶[i] = g✶[j]
    end
    # branch2_g✶, branch2_f, branch1_g✶, branch1_f
    # branch1_g✶, branch1_f, branch2_g✶, branch2_f
    # mid1 
    if branch1_f[2][2] > branch1_f[1][2]
        branch2_g✶, branch2_f, branch1_g✶, branch1_f
    else
        branch1_g✶, branch1_f, branch2_g✶, branch2_f
    end
end

function isoredshift!(data::Vector{<:RedshiftData}, g✶, which = :lower)
    vals = if which == :lower
        map(i -> interpolate(i, g✶, false), data)
    elseif which == :upper
        map(i -> interpolate(i, g✶, true), data)
    else
        error("bad symbol")
    end
    vals = reduce(hcat, vals)
    vals[1, :], vals[2, :]
end

function plot_branches_no_text(ax, rf; kwargs...)
    a = copy(rf.α)
    b = copy(rf.β)
    a[@.(!rf.visible)] .= NaN
    b[@.(!rf.visible)] .= NaN
    lines!(ax, a, b; kwargs...)
end

function plot_branches(
    ax,
    rf::RedshiftData,
    annotate,
    add_label = false;
    color = :black,
    kwargs...,
)
    plot_branches_no_text(ax, rf; color = color, kwargs...)
    x2 = rf.upper.α_interp.u
    y2 = rf.upper.β_interp.u
    r_str = Printf.@sprintf("%0.0f", rf.r)
    text = if add_label
        L"r_\text{em} = %$(r_str) \, r_\text{g}"
    else
        L"%$(r_str)"
    end
    if annotate
        text!(
            ax,
            x2[1] + 0.2,
            y2[1] - 0.5,
            text = text,
            color = color,
            fontsize = SMALL_FONT_SIZE,
        )
    end
end

function calculate_redshift_data(m, x, d, r; kwargs...)
    rshift = ConstPointFunctions.redshift(m, x)
    α, β = impact_parameters_for_radius(m, x, d, r; N = 300, kwargs...)
    vs = map_impact_parameters.(m, (x,), α, β)
    points = tracegeodesics(
        m,
        fill(x, size(vs)),
        vs,
        d isa AbstractThickAccretionDisc ?
        Gradus.DatumPlane(Gradus.cross_section(d, r)) : d,
        x[2] * 2,
        ensemble = Gradus.EnsembleEndpointThreads(),
        save_on = false,
    )
    g = [rshift(m, gp, 0.0) for gp in points]
    points2 = tracegeodesics(
        m,
        fill(x, size(vs)),
        vs,
        d,
        x[2] * 2,
        ensemble = Gradus.EnsembleEndpointThreads(),
        save_on = false,
    )
    is_visible(gp1, gp2) = abs(gp1.x[2] * sin(gp1.x[3]) - gp2.x[2] * sin(gp2.x[3])) < 1e-2
    I = @. is_visible(points2, points)
    RedshiftData(α, β, g, r, I)
end

function plot_parameterization!(ax, data, radial, minor_radial, g✶s, palette)
    for (i, d) in enumerate(minor_radial)
        plot_branches_no_text(
            ax,
            d,
            color = :lightgrey,
            linewidth = 1.1,
            # alpha = 0.5,
        )
    end

    # N = 10
    # offset = trunc(Int, length(radii) / N)
    # I = trunc.(Int, range(offset, length(radii) - offset, 5))
    # for d in data[I]
    for (i, d) in enumerate(radial)
        plot_branches(
            ax,
            d,
            true,
            i == lastindex(radial),
            color = popfirst!(palette),
            linewidth = 1.5,
        )
    end

    for g in g✶s
        color = popfirst!(palette)
        # don't use yellow in this plot since it's tricky to see
        if color.b == 0.25882354f0
            color = popfirst!(palette)
        end
        α1, β1 = isoredshift!(data, g, :upper)
        α2, β2 = isoredshift!(data, g, :lower)
        lines!(ax, α1, β1, color = color, linewidth = 2.0)
        lines!(ax, α2, β2, color = color, linewidth = 2.0)
        text = if g == g✶s[end-1]
            L"g^\ast=%$(g)"
        else
            L"%$(g)"
        end
        text!(
            ax,
            α1[end] - (g == g✶s[end-1] ? 0.8 : 0.0),
            β1[end] - 0.2,
            text = text,
            align = (:center, :top),
            color = color,
            fontsize = SMALL_FONT_SIZE,
        )
    end

    plot_branches(ax, data[1], false, color = :black, linewidth = 3)
end

function calculate_disc_parameterization(m, x, d; kwargs...)
    radii = sort(vcat(range(Gradus.isco(m) + 1e-2, 5, 20), range(5, 15.0, 20)))
    # radii = range(Gradus.isco(m) + 1e-2, 5, 20)
    data = @time Gradus._threaded_map(radii) do r
        calculate_redshift_data(m, x, d, r; kwargs...)
    end
    selected_radii = [2.0, 5.0, 8.0, 11.0, 14.0]
    radial = @time Gradus._threaded_map(selected_radii) do r
        calculate_redshift_data(m, x, d, r; kwargs...)
    end
    minor_radii = filter(>=(Gradus.isco(m)), range(-1, 5, step = 0.5))
    minor_radial = @time Gradus._threaded_map(minor_radii) do r
        calculate_redshift_data(m, x, d, r; kwargs...)
    end

    F(collection) = filter(!isnothing, collection)
    F(data), F(radial), F(minor_radial)
end

m = KerrMetric(1.0, 0.998)
x = SVector(0.0, 10_000.0, deg2rad(75), 0.0)
d = ThinDisc(0.0, 1000.0)

g✶s = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 0.95]

data1, radial1, minor_radial1 = calculate_disc_parameterization(m, x, d)

d2 = ShakuraSunyaev(m)
data2, radial2, minor_radial2 = calculate_disc_parameterization(m, x, d2, β₀ = 1)

begin
    fig = Figure(resolution = (1200, 300))
    ax1 = Axis(
        fig[1, 1],
        aspect = DataAspect(),
        ylabel = L"\beta",
        xlabel = L"\alpha",
        xticks = [-15, -10, -5, 0, 5, 10, 15],
        topspinevisible = false,
        leftspinevisible = false,
        rightspinevisible = false,
        bottomspinevisible = false,
    )
    ax2 = Axis(
        fig[1, 2],
        aspect = DataAspect(),
        ylabel = L"\beta",
        xlabel = L"\alpha",
        xticks = [-15, -10, -5, 0, 5, 10, 15],
        topspinevisible = false,
        leftspinevisible = false,
        rightspinevisible = false,
        bottomspinevisible = false,
    )
    hidedecorations!(ax1)
    hidedecorations!(ax2)

    palette = _default_palette()
    plot_parameterization!(ax1, data1, radial1, minor_radial1, g✶s, palette)
    palette = _default_palette()
    plot_parameterization!(ax2, data2, radial2, minor_radial2, g✶s, palette)

    linkyaxes!(ax1, ax2)
    xlims!(ax1, -19, 22)
    xlims!(ax2, -19, 22)
    ylims!(ax2, -6, nothing)

    text!(ax1, -13, 10, text = "a", fontsize = LARGE_FONT_SIZE, font = :bold)
    text!(ax2, -13, 10, text = "b", fontsize = LARGE_FONT_SIZE, font = :bold)

    @savefigure(fig)
    fig
end
