#=
figure_3.jl

Draws Figure 3 (Section 5), figures/figure_3.pdf, for Proposition 2 and
the inner boundary (Proposition 3). Along the ray -μ/σ² = λ in the parent
parameters of the normal, the truncation to I is ∝ exp(-λx - x²/(2σ²)), and
as σ → ∞ it tends to the truncated exponential with rate λ.
(a) Schematic of the rays μ = -λσ² in the parent parameters, for
    λ = -2, 0.5, 3, as arrows in the direction of increasing σ; legend here.
(b) I = [0,1]: the (mean, sd) trajectories along these rays, each ending at
    its point of Γ_I, with the inner boundary Γ_I^in (dashed; the μ = 0
    trajectory and its mirror image) and the inner region shaded.
(c) I = [0,∞): the trajectories along the rays λ = 0.5 (outer, ending on
    the line sd = mean) and λ = -2 (inner, moving away from it), with the
    inner boundary sd = sqrt(π/2-1) mean (dashed).
Layout: (a) and (b) side by side, (c) centred below them. Dashed lines mark
the borders between inner and outer truncations in every panel.
Markers of the same shape in all panels are at the same σ².
=#

using Plots, LaTeXStrings, QuadGK
gr()

const BLUE = RGB(0.16, 0.24, 0.55)
const COLS = (RGB(0.90, 0.50, 0.05), RGB(0.10, 0.55, 0.20), RGB(0.70, 0.10, 0.45))
const LAMS = (-2.0, 0.5, 3.0)                     # rays in panel (b)
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
# The vertical lines μ = 0 (the ray λ = 0) and μ = 1 bound the inner
# truncations on [0,1]; they carry the only two ticks.
pa = plot(; xlims = (-7, 5), ylims = (0, 2.1), framestyle = :axes, grid = false,
          xticks = ([0, 1], ["0", "1"]), yticks = false, legend = :bottomleft,
          background_color_legend = :white, foreground_color_legend = :gray70,
          xlabel = L"\mu", ylabel = L"\sigma^2", title = "(a)", FONT...)
plot!(pa, [0, 0], [0, 2.1]; color = :gray70, lw = 0.8, ls = :dash, label = "")
plot!(pa, [1, 1], [0, 2.1]; color = :gray70, lw = 0.8, ls = :dash, label = "")
for (λ, col) in zip(LAMS, COLS)
    plot!(pa, [0, -2λ], [0, 2]; color = col, lw = 1.6, label = L"\lambda={%$(isinteger(λ) ? Int(λ) : λ)}",
          arrow = arrow(:closed, :head, 0.4, 0.3))
    scatter!(pa, [-λ * s for s in S2], collect(S2); color = col, label = "",
             markershape = collect(SHAPES), markersize = 3.5, markerstrokewidth = 0)
end

# The inner boundary on [0,1]: the μ = 0 trajectory (the ray λ = 0) and its
# mirror image under x ↦ 1 - x.
inn = vcat([(0.0, 0.0)], [moments(0.0, s) for s in exp10.(range(-5, 4; length = 400))])
inn = vcat(inn, reverse([(1 - m, s) for (m, s) in inn]))

# (b) trajectories in the (mean, sd) plane, with Γ_I and Γ_I^in.
λs = vcat(-reverse(exp10.(range(-3, log10(300); length = 300))), 0.0,
          exp10.(range(-3, log10(300); length = 300)))
gam = [moments(λ, Inf) for λ in λs]
pb = plot(first.(gam), last.(gam); color = BLUE, lw = 2.2, label = "",
          xlims = (0, 1), ylims = (0, 0.32), framestyle = :box, legend = false,
          xlabel = L"\bar\mu", ylabel = L"\bar\sigma", title = "(b)", gridalpha = 0.08,
          background_color_legend = :white, foreground_color_legend = :gray70, FONT...)
plot!(pb, first.(inn), last.(inn); fillrange = 0, fillcolor = :gray88, fillalpha = 1, color = :gray88,
      lw = 0, label = "")
plot!(pb, first.(inn), last.(inn); color = :gray25, lw = 1.2, ls = :dash, label = "")
plot!(pb, first.(gam), last.(gam); color = BLUE, lw = 2.2, label = "")
for (λ, col) in zip(LAMS, COLS)
    # As σ → 0 the truncation concentrates at 0, so each trajectory starts at
    # the origin.
    traj = vcat([(0.0, 0.0)], [moments(λ, s) for s in exp10.(range(-5, 4; length = 500))])
    plot!(pb, first.(traj), last.(traj); color = col, lw = 1.6, label = "")
    mk = [moments(λ, s) for s in S2]
    scatter!(pb, first.(mk), last.(mk); color = col, markershape = collect(SHAPES),
             markersize = 3.5, markerstrokewidth = 0, label = "")
end
# (c) the half-line [0,∞): the rays λ = 0.5 (outer, limiting mean 2) and
# λ = -2 (inner, μ = 2σ² > 0, coefficient of variation decreasing to 0). The
# ray λ = 3 is a scaled copy of λ = 0.5 and is omitted to keep the panel clear.
function moments_half(λ, σ2)
    c = λ < 0 ? λ^2 * σ2 / 2 : 0.0                 # log-density maximum, against overflow
    g = x -> exp(-λ * x - x^2 / (2σ2) - c)
    Z = quadgk(g, 0, Inf)[1]
    m = quadgk(x -> x * g(x), 0, Inf)[1] / Z
    return m, sqrt(quadgk(x -> (x - m)^2 * g(x), 0, Inf)[1] / Z)
end
const CIN = sqrt(π / 2 - 1)
# Drawn like Figure 1 (b): the half-line region is an unbounded cone, so no
# frame or ticks; the shaded inner cone ends before the right edge.
pc = plot([0, 3.2, 3.2, 0], [0, 3.2 * CIN, 0, 0]; seriestype = :shape, fillcolor = :gray88,
          linecolor = :gray88, xlims = (0, 3.4), ylims = (0, 2.55), legend = false,
          framestyle = :axes, grid = false, xticks = false, yticks = false,
          xlabel = L"\bar\mu", ylabel = L"\bar\sigma", title = "(c)",
          right_margin = 15Plots.mm, FONT...)
plot!(pc, [0, 3.2], [0, 3.2 * CIN]; color = :gray25, lw = 1.2, ls = :dash, label = "")
plot!(pc, [0, 2.45], [0, 2.45]; color = BLUE, lw = 2.2, label = "",
      arrow = arrow(:closed, :head, 0.4, 0.3))
annotate!(pc, 1.0, 1.35, text(L"\bar\sigma=\bar\mu", 9, BLUE, :right, "Computer Modern"))
for (λ, col, s2max) in ((0.5, COLS[2], 1e5), (-2.0, COLS[1], 1.57))
    traj = vcat([(0.0, 0.0)], [moments_half(λ, s) for s in exp10.(range(-5, log10(s2max); length = 500))])
    plot!(pc, first.(traj), last.(traj); color = col, lw = 1.6, label = "")
    mk = [moments_half(λ, s) for s in S2]
    scatter!(pc, first.(mk), last.(mk); color = col, markershape = collect(SHAPES),
             markersize = 3.5, markerstrokewidth = 0, label = "")
end
# (a) and (b) side by side, (c) centred below them. The right margin of (c)
# balances its axis labels on the left, so that its frame sits in the middle.
fig = plot(pa, pb, pc;
           layout = @layout([a b; _ c{0.58w} _]), size = (620, 480),
           left_margin = 5Plots.mm, bottom_margin = 4Plots.mm)
savefig(fig, joinpath(@__DIR__, "..", "figures", "figure_3.pdf"))
println("wrote figures/figure_3.pdf")
