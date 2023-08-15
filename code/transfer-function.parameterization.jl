using Makie, CairoMakie, LaTeXStrings, Printf
using Gradus

include("common.jl")

struct RedshiftBranch{K}
    α_interp::K
    β_interp::K
end

function RedshiftBranch(g✶, ps)
    I = sortperm(g✶)
    α = first.(ps)[I]
    β = last.(ps)[I]
    t1 = Gradus.DataInterpolations.LinearInterpolation(α, g✶[I])
    t2 = Gradus.DataInterpolations.LinearInterpolation(β, g✶[I])
    RedshiftBranch(t1, t2)
end

struct RedshiftData{T,K}
    upper::RedshiftBranch{K}
    lower::RedshiftBranch{K}
    gmin::T
    gmax::T
    r::T
end

function RedshiftData(α, β, g, r)
    gmin, gmax = extrema(g)
    g✶ = Gradus.g_to_g✶.(g, gmin, gmax)

    points = [SVector(a, b) for (a, b) in zip(α, β)]
    gup, pup, gdown, pdown = splitbranches(g✶, points)

    RedshiftData(RedshiftBranch(gup, pup), RedshiftBranch(gdown, pdown), gmin, gmax, r)
end

function interpolate(rb::RedshiftBranch, g✶)
    SVector(rb.α_interp(g✶), rb.β_interp(g✶))
end

function splitbranches(g✶::Vector{T}, f::Vector{V}) where {T,V}
    _, imin = findmin(g✶)
    _, imax = findmax(g✶)
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

    if branch1_f[2][2] > branch1_f[1][2]
        branch2_g✶, branch2_f, branch1_g✶, branch1_f
    else
        branch1_g✶, branch1_f, branch2_g✶, branch2_f
    end
end

function isoredshift!(data::Vector{<:RedshiftData}, g✶, which = :lower)
    vals = if which == :lower
        map(i -> interpolate(i.lower, g✶), data)
    elseif which == :upper
        map(i -> interpolate(i.upper, g✶), data)
    else
        error("bad symbol")
    end
    vals = reduce(hcat, vals)
    vals[1, :], vals[2, :]
end

function plot_branches(
    ax,
    rf::RedshiftData,
    annotate,
    add_label = false;
    color = :black,
    kwargs...,
)
    x1 = rf.lower.α_interp.u
    y1 = rf.lower.β_interp.u
    lines!(ax, x1, y1; color = color, kwargs...)
    x2 = rf.upper.α_interp.u
    y2 = rf.upper.β_interp.u
    lines!(ax, x2, y2; color = color, kwargs...)
    r_str = Printf.@sprintf("%0.0f", rf.r)
    text = if add_label
        L"r_\text{em} = %$(r_str)"
    else
        L"%$(r_str)"
    end
    if annotate
        text!(ax, x2[1] + 0.2, y2[1], text = text, color = color)
    end
end

function calculate_redshift_data(m, x, d, r)
    rshift = ConstPointFunctions.redshift(m, x)
    α, β = impact_parameters_for_radius(m, x, d, r, N = 500)
    vs = map_impact_parameters.(m, (x,), α, β)
    points = tracegeodesics(
        m,
        fill(x, size(vs)),
        vs,
        d,
        x[2] * 2,
        ensemble = Gradus.EnsembleEndpointThreads(),
        save_on = false,
    )
    g = [rshift(m, gp, 0.0) for gp in points]
    RedshiftData(α, β, g, r)
end


m = KerrMetric(1.0, 0.998)
x = SVector(0.0, 10_000.0, deg2rad(75), 0.0)
d = GeometricThinDisc(0.0, 1000.0, π / 2)

radii = range(Gradus.isco(m) + 1e-2, 16, 34)
data = @time Gradus._threaded_map(radii) do r
    calculate_redshift_data(m, x, d, r)
end

selected_radii = [2.0, 5.0, 9.0, 12.0, 15.0]
radial = @time Gradus._threaded_map(selected_radii) do r
    calculate_redshift_data(m, x, d, r)
end

g✶s = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 0.95]

begin
    palette = Iterators.Stateful(Iterators.Cycle(Makie.wong_colors()))
    fig = Figure(resolution = (600, 280))
    ax = Axis(
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
    hidedecorations!(ax)

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
        )
    end

    # xlims!(ax, nothing, 0.0)
    ylims!(ax, -6, nothing)
    xlims!(ax, -19, 21)

    plot_branches(ax, data[1], false, color = :black, linewidth = 3)

    fig
    @savefigure(fig)
end
