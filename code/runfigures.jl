macro buildfigure(path)
    quote
        @info "Building figure: " * $(path)
        @time include($path)
    end |> esc
end

# on my laptop it takes about 20 minutes to build all of the figures 
# so set this script running, go for a coffee, and then come back
# and compile the LaTeX document :)

@buildfigure("circular-orbits.E-Lz.jl")
@buildfigure("deflection.iyer-hansen.jl")
@buildfigure("emissivity.coronal-traces.jl")
@buildfigure("emissivity.point-source.jl")
@buildfigure("lineprofiles.comparison.jl")
@buildfigure("lineprofiles.ssd.jl")
@buildfigure("radiative-transfer.gold.jl")
@buildfigure("reverberation.lag-energy.jl")
@buildfigure("reverberation.thin-disc.jl")
@buildfigure("skycoords.jl")
@buildfigure("stability.conservation.jl")
@buildfigure("transfer-function.parameterization.jl")
@buildfigure("transfer-functions.2d.jl")
@buildfigure("transfer-functions.plots.jl")
