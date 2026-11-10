#=
figure_3.jl

Draws Figure 3, the illustration of Proposition 2, figures/figure_3.pdf, for
I = [0,1]. Along the ray -μ/σ² = λ in the parent parameters of the normal,
the truncation to [0,1] is ∝ exp(-λx - x²/(2σ²)), and as σ → ∞ it tends to
the truncated exponential with rate λ.
(a) Schematic of the rays μ = -λσ² in the parent parameters, for
    λ = -2, 3, 8, as arrows in the direction of increasing σ; legend here.
(b) The (mean, sd) trajectories along these rays, each ending at its point
    of Γ_I (open circle). Markers of the same shape in both panels are at
    the same σ².
=#

using Plots, LaTeXStrings, QuadGK
gr()

const BLUE = RGB(0.16, 0.24, 0.55)
const COLS = (RGB(0.90, 0.50, 0.05), RGB(0.10, 0.55, 0.20), RGB(0.70, 0.10, 0.45))
const LAMS = (-2.0, 3.0, 8.0)                     # rays in panel (b)
const S2 = (0.1, 0.5, 1.5)                        # σ² values marked in both panels
const SHAPES = (:circle, :square, :utriangle)
const FONT = (fontfamily = "Computer Modern", titlefontsize = 10, guidefontsize = 10,
              tickfontsize = 8, legendfontsize = 8)

lw(λ, σ2) = x -> -λ * x - x^2 / (2σ2)             # log of the unnormalised density
function moments(λ, σ2)
    g = x -> exp(lw(λ, σ2)(x) - max(0, -λ))       # shifted against overflow
    Z = quadgk(g, 0, 1)[1]
    m = quadgk(x -> x * g(x), 0, 1)[1] / Z
    return m, sqrt(quadgk(x -> (x - m)^2 * g(x), 0, 1)[1] / Z)
end
# (a) schematic of the rays in the parent parameters, drawn as arrows in the
# direction of increasing σ; the legend for both panels is here.
pa = plot(; xlims = (-17, 5), ylims = (0, 2.1), framestyle = :axes, grid = false,
          xticks = false, yticks = false, legend = :bottomleft,
          background_color_legend = :white, foreground_color_legend = :gray70,
          xlabel = L"\mu", ylabel = L"\sigma^2", title = "(a)", FONT...)
plot!(pa, [0, 0], [0, 2.1]; color = :gray70, lw = 0.8, label = "")
annotate!(pa, 0.25, 2.0, text(L"\mu=0", 8, :left, :gray40, "Computer Modern"))
for (λ, col) in zip(LAMS, COLS)
    plot!(pa, [0, -2λ], [0, 2]; color = col, lw = 1.6, label = L"\lambda=%$(Int(λ))",
          arrow = arrow(:closed, :head, 0.4, 0.3))
    scatter!(pa, [-λ * s for s in S2], collect(S2); color = col, label = "",
             markershape = collect(SHAPES), markersize = 3.5, markerstrokewidth = 0)
end

# (b) trajectories in the (mean, sd) plane, with Γ_I.
λs = vcat(-reverse(exp10.(range(-3, log10(300); length = 300))), 0.0,
          exp10.(range(-3, log10(300); length = 300)))
gam = [moments(λ, Inf) for λ in λs]
pb = plot(first.(gam), last.(gam); color = BLUE, lw = 2.2, label = "",
          xlims = (0, 0.7), ylims = (0, 0.32), framestyle = :box, legend = false,
          xlabel = L"\bar\mu", ylabel = L"\bar\sigma", title = "(b)", gridalpha = 0.08,
          background_color_legend = :white, foreground_color_legend = :gray70, FONT...)
for (λ, col) in zip(LAMS, COLS)
    traj = [moments(λ, s) for s in exp10.(range(-3, 4; length = 400))]
    plot!(pb, first.(traj), last.(traj); color = col, lw = 1.6, label = "")
    mk = [moments(λ, s) for s in S2]
    scatter!(pb, first.(mk), last.(mk); color = col, markershape = collect(SHAPES),
             markersize = 3.5, markerstrokewidth = 0, label = "")
    m∞, s∞ = moments(λ, Inf)
    scatter!(pb, [m∞], [s∞]; color = :white, markerstrokecolor = col,
             markerstrokewidth = 1.5, markersize = 4.5, label = "")
end
fig = plot(pa, pb; layout = (1, 2), size = (620, 240), left_margin = 5Plots.mm,
           bottom_margin = 3Plots.mm)
savefig(fig, joinpath(@__DIR__, "..", "figures", "figure_3.pdf"))
println("wrote figures/figure_3.pdf")
