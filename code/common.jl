using Makie, CairoMakie, LaTeXStrings

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
