using Makie, CairoMakie, LaTeXStrings

# this function is absolutely terrible but it does the job of drawing an eye
function draw_observer_eye!(ax, x0, y0, scale; flip = false, linewidth = 4.0, rot = 0)
    rotmat = [cos(rot) -sin(rot) ; sin(rot) cos(rot)]

    x_xfm(x) = @. scale * x * (flip ? -1 : 1) + x0
    y_xfm(y) = @. scale * y + y0

    function xfm(x, y)
        rotated = map(eachindex(x)) do i
            rotmat * [x[i], y[i]]
        end
        xr = first.(rotated)
        yr = last.(rotated)
        x_xfm(xr), y_xfm(yr)
    end

    ls = collect(range(0, 2, 20))
    x1 = ls
    y1 = @. exp(0.3 * ls) - 1 

    kwargs = (; color = :black, linewidth = linewidth)
    lines!(ax, xfm(x1, y1)...; kwargs...)
    lines!(ax, xfm(x1[1:end-1], -y1[1:end-1])...; kwargs...)

    i1 = 2
    θ1 = atan(y1[end-i1], x1[end-1])
    θ2 = atan(-y1[end-i1], x1[end-1])

    theta = collect(range(θ1, θ2, 20))

    x2 = @. x1[end-i1] * cos(theta)
    y2 = @. x1[end-i1] * sin(theta)

    lines!(ax, xfm(x2, y2)...; kwargs...)

    i2 = 1
    phi1 = atan(y2[i2], x2[i2])
    phi2 = atan(y2[end-i2], x2[end-i2])

    phi = collect(range(phi1, phi2, 20))
    x3 = @. -1.7 * cos(phi) + (1.7 + x1[end-i1] - 0.15) 
    y3 = @. 1 * sin(phi) 

    lines!(ax, xfm(x3, y3)...; kwargs...)
end

FIGURE_DIR = joinpath(@__DIR__(), "..", "latex", "figures")

macro savefigure(fig)
    quote
        _save_figure(fig, splitpath($(string(__source__.file)))[end])
    end
end

function _default_palette()
    Iterators.Stateful(Iterators.Cycle(Makie.wong_colors()))
end

function _save_figure(fig, filename)
    if filename[end-2:end] == ".jl"
        filename = filename[1:end-3] * ".pdf"
    end
    path = joinpath(FIGURE_DIR, filename)
    @info "Saving figure to : $(path)"
    Makie.save(joinpath(FIGURE_DIR, filename), fig)
end