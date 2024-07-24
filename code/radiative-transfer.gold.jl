using CairoMakie, Makie, LaTeXStrings, Printf
using Gradus

include("common.jl")

# Eq (5)
function number_density(r, cosθ, h, n₀)
    z = h * cosθ
    n₀ * exp(-0.5 * ((r / 10)^2 + z^2))
end

# Eq (6)
function angular_momentum(r, sinθ, l₀, q = 0.5)
    R = r * sinθ
    (l₀ / (1 + R)) * R^(1 + q)
end

# Eq (8)
function ū(g, l)
    √(-(inv(g[1] - 2g[5] * l + g[4] * l^2)))
end

# Eq (7)
function four_velocity(m::AbstractMetric, r, θ, l)
    ginv = Gradus.inverse_metric_components(m, SVector(r, θ))
    M = Gradus._symmetric_matrix(ginv)
    A = ū(ginv, l)
    M * SVector(-A, 0, 0, A * l)
end

# Eq (9)
function fluid_emissivity(C, n, ν, α; νp = 1)
    C * n * (ν / νp)^(-α)
end

# Eq (11)
function fluid_absorptivity(C, n, ν, α, A; β = 2.5, νp = 1)
    A * C * n * (ν / νp)^(-(β + α))
end

function Sν(ν, νp, A, β)
    inv(A) * (ν / νp)^β
end


struct AnalyticDiscTest{T} <: AbstractAccretionDisc{T}
    A::T
    α::T
    h::T
    l₀::T
end
AnalyticDiscTest(; A = 0.0, α = -3.0, h = 0.0, l₀ = 0.0) = AnalyticDiscTest(A, α, h, l₀)
Gradus.is_finite_disc(::Type{<:AnalyticDiscTest}) = false

function Gradus.fluid_velocity(m::AbstractMetric, d::AnalyticDiscTest, x, r_isco, λ)
    sinθ = sin(x[3])
    l = angular_momentum(x[2], sinθ, d.l₀)
    four_velocity(m, x[2], x[3], l)
end

function Gradus.fluid_absorption_emission(::AbstractMetric, d::AnalyticDiscTest, x, ν, u)
    # i really don't know what to set for these
    # "physical values" so i am just setting C to make sure the 
    # flux is approximately correct for one of them
    C = 6.2104e-6 
    n₀ = 1
    νp = 1
    cosθ = cos(x[3])
    n = number_density(x[2], cosθ, d.h, n₀)

    jν = fluid_emissivity(C, n, ν, d.α; νp = νp)
    αν = fluid_absorptivity(C, n, ν, d.α, d.A; νp = νp)

    αν, jν
end

function do_trace(m, x, disc)
    @time rendergeodesics(
        m,
        x,
        disc,
        2x[2],
        verbose = true,
        αlims = (-15, 15),
        βlims = (-15, 15),
        image_width = 128,
        image_height = 128,
        trace = Gradus.TraceRadiativeTransfer(I₀ = 0.0),
        pf = PointFunction((m, gp, t) -> gp.aux[1]),
        chart = Gradus.chart_for_metric(m, 2x[2]),
    )
end

M = 1.0
m = KerrMetric(M, 0.9)
x = SVector(0.0, 1000.0, deg2rad(60), 0.0)

# same with A, since I don't really know how to convert C I don't know how to convert A, so I estimate
A_fudge = 0.289

test1 = AnalyticDiscTest(A = 0.0, α = -3.0, h = 0.0, l₀ = 0.0)
test2 = AnalyticDiscTest(A = 0.0, α = -2.0, h = 0.0, l₀ = 1.0)
m_test2 = KerrMetric(1.0, 0.0)
test3 = AnalyticDiscTest(A = 0.0, α = 0.0, h = 10.0 / 3, l₀ = 1.0)
# similarly here, not sure what A should be set to, but there must
# be an order magnitude difference between these two
test4 = AnalyticDiscTest(A = A_fudge * 1e5, α = 0.0, h = 10.0 / 3, l₀ = 1.0)
test5 = AnalyticDiscTest(A = A_fudge * 1e6, α = 0.0, h = 100.0 / 3, l₀ = 1.0)

a1, b1, im1 = do_trace(m, x, test1)
a2, b2, im2 = do_trace(m_test2, x, test2)
a3, b3, im3 = do_trace(m, x, test3)
a4, b4, im4 = do_trace(m, x, test4)
a5, b5, im5 = do_trace(m, x, test5)

function plot_heatmap!(ax, x, y, im; kwargs...)
    S = sum(filter(!isnan, im))
    hm = heatmap!(ax, x, y, im' ./ S, colormap = :cubehelix; kwargs...)
    s = Printf.@sprintf("%1.4f", S) #/ 1e5)
    text!(
        ax,
        (-14.5, 10),
        text = L"S_\text{tot} = %$(s)\, \text{Jy}",
        fontsize = 24,
        color = :white,
    )
    hm
end

begin
    fig = Figure(resolution = (1100, 300))

    ga = fig[1, 1] = GridLayout()

    ax1 = Axis(ga[1, 1], aspect = DataAspect(), ylabel = L"\beta", xlabel = L"\alpha")
    hm1 = plot_heatmap!(ax1, a1, b1, im1)
    crange = (0.0, 2e-4)
    Label(ga[1, 1, Top()], "Test 1", padding = (0, 0, 0, 0), font = :bold)

    ax2 = Axis(ga[1, 2], aspect = DataAspect(), xlabel = L"\alpha")
    hm2 = plot_heatmap!(ax2, a2, b2, im2, colorrange = crange)
    Label(ga[1, 2, Top()], "Test 2", padding = (0, 0, 0, 0), font = :bold)

    ax3 = Axis(ga[1, 3], aspect = DataAspect(), xlabel = L"\alpha")
    hm3 = plot_heatmap!(ax3, a3, b3, im3, colorrange = crange)
    Label(ga[1, 3, Top()], "Test 3", padding = (0, 0, 0, 0), font = :bold)

    ax4 = Axis(ga[1, 4], aspect = DataAspect(), xlabel = L"\alpha")
    hm4 = plot_heatmap!(ax4, a4, b4, im4, colorrange = crange)
    Label(ga[1, 4, Top()], "Test 4", padding = (0, 0, 0, 0), font = :bold)

    ax5 = Axis(ga[1, 5], aspect = DataAspect(), xlabel = L"\alpha")
    hm5 = plot_heatmap!(ax5, a5, b5, im5, colorrange = crange)
    Label(ga[1, 5, Top()], "Test 5", padding = (0, 0, 0, 0), font = :bold)

    dp1(x) = Printf.@sprintf("%.1f", x)
    cbar = Colorbar(
        ga[1, 6],
        limits = crange,
        label = L"S / S_\text{tot}",
        labelsize = 24,
        ticks = [0, 0.5e-4, 1.0e-4, 1.5e-4, 2e-4],
        tickformat = values -> [L"%$(dp1(x * 1e4)) \times 10^{-4}" for x in values],
        tellheight = false,
        colormap = :cubehelix,
        height = 172,
    )

    linkyaxes!(ax1, ax2, ax3, ax4, ax5)
    hideydecorations!(ax2)
    hideydecorations!(ax3)
    hideydecorations!(ax4)
    hideydecorations!(ax5)
    colgap!(ga, 8)

    resize_to_layout!(fig)
    @savefigure(fig)
    fig
end
