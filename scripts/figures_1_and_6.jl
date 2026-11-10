#=
figures_1_and_6.jl

Draws Figure 1 (figures/figure_1.pdf) and Figure 6 (figures/figure_6.pdf).

The boundary Γ_I^sd for I = [0,1] is traced by the truncated exponential
∝ exp(-λx) on [0,1] as λ ranges over ℝ. The conjectured Student-t_ν domes
are traced in the same way by the tail-tilt family
∝ (1 + λx/(ν+1))^{-(ν+1)}, for λ ≥ 0 and mirrored about 1/2.
=#

using Plots, LaTeXStrings, QuadGK
gr()

const OUT = joinpath(@__DIR__, "..", "figures")
const BLUE = RGB(0.16, 0.24, 0.55)
const FONT = (fontfamily = "Computer Modern", titlefontsize = 10, guidefontsize = 10,
              tickfontsize = 8, legendfontsize = 7)

# Mean and sd of the density ∝ g(x) on [0,1], for g decreasing (mean ≤ 1/2).
function moments(g)
    Z = quadgk(g, 0, 1)[1]
    m = quadgk(x -> x * g(x), 0, 1)[1] / Z
    return m, sqrt(quadgk(x -> (x - m)^2 * g(x), 0, 1)[1] / Z)
end

# Trace a symmetric arch from a one-parameter decreasing family g_λ, λ ≥ 0.
function arch(g)
    λs = vcat(0.0, exp10.(range(-3, 5; length = 400)))
    pts = [moments(g(λ)) for λ in λs]
    left = vcat([(0.0, 0.0)], reverse(pts))      # from the corner to the top
    right = [(1 - m, s) for (m, s) in reverse(left)]
    pts = vcat(left, right[2:end])
    return first.(pts), last.(pts)
end

exparch() = arch(λ -> (x -> exp(-λ * x)))
tarch(ν) = arch(λ -> (x -> (1 + λ * x / (ν + 1))^(-(ν + 1))))

# Figure 1.
mg, sg = exparch()
p1 = plot(mg, sg; fillrange = 0, fillcolor = :gray85, linecolor = BLUE, lw = 2.2,
          xlims = (-0.01, 1.01), ylims = (0, 0.52), framestyle = :box, legend = false,
          xlabel = L"\bar\mu", ylabel = L"\bar\sigma", title = "(a)", gridalpha = 0.08,
          FONT...)
scatter!(p1, [0.5], [1 / sqrt(12)]; color = :black, markersize = 4, markerstrokewidth = 0)
# The Bhatia-Davis bound sd = sqrt(m(1-m)), valid for every distribution on [0,1].
mb = range(0, 1; length = 401)
plot!(p1, mb, sqrt.(mb .* (1 .- mb)); color = :gray45, lw = 1.0)
# The grey cone ends before the arrowhead.
p2 = plot([0, 9.5, 9.5, 0], [0, 9.5, 0, 0]; seriestype = :shape, fillcolor = :gray85,
          linecolor = :gray85, xlims = (0, 10.5), ylims = (0, 10.5), legend = false,
          framestyle = :axes, grid = false, xticks = false, yticks = false,
          xlabel = L"\bar\mu", ylabel = L"\bar\sigma", title = "(b)", FONT...)
plot!(p2, [0, 10.2], [0, 10.2]; color = BLUE, lw = 2.2, arrow = arrow(:closed, :head, 0.4, 0.3))
annotate!(p2, 5.3, 6.6, text(L"\bar\sigma=\bar\mu", 10, BLUE, :right, "Computer Modern"))
fig1 = plot(p1, p2; layout = (1, 2), size = (620, 240), left_margin = 2Plots.mm,
            bottom_margin = 3Plots.mm)
savefig(fig1, joinpath(OUT, "figure_1.pdf"))

# Figure 4.
# Solid lines in distinct colours (orange, red, purple) for nu = 3, 5, 10.
styles = ((3, RGB(0.90, 0.50, 0.05), :solid), (5, RGB(0.80, 0.10, 0.15), :solid),
          (10, RGB(0.50, 0.20, 0.65), :solid))
q1 = plot(; xlims = (-0.01, 1.01), ylims = (0, 0.32), framestyle = :box,
          xlabel = L"\bar\mu", ylabel = L"\bar\sigma", title = "(a)",
          legend = :bottom, background_color_legend = :white,
          foreground_color_legend = :gray70, gridalpha = 0.08, FONT...)
plot!(q1, mg, sg; fillrange = 0, fillcolor = :gray85, fillalpha = 0.8, linecolor = BLUE,
      lw = 2.2, label = "Log-concave")
for (ν, col, ls) in styles
    m, s = tarch(ν)
    plot!(q1, m, s; color = col, ls = ls, lw = 1.8, label = "Student-\$t_{$ν}\$")
end
q2 = plot(; xlims = (0, 10), ylims = (0, 10), framestyle = :box,
          xlabel = L"\bar\mu", ylabel = L"\bar\sigma", title = "(b)",
          legend = :bottomright, background_color_legend = :white,
          foreground_color_legend = :gray70, gridalpha = 0.08, FONT...)
plot!(q2, [0, 10], [0, 10]; color = BLUE, lw = 2.2,
      label = "Log-concave: \$\\bar\\sigma\\leq\\bar\\mu\$")
for (ν, col, ls) in styles
    c = sqrt(ν / (ν - 2))
    r = ν == 3 ? "3" : "$ν/$(ν-2)"
    lab = L"t_{%$ν}\!:\ \bar\sigma\leq\sqrt{%$r}\,\bar\mu"
    plot!(q2, [0, 10], [0, 10c]; color = col, ls = ls, lw = 1.8, label = lab)
end
fig4 = plot(q1, q2; layout = (1, 2), size = (620, 250), left_margin = 2Plots.mm,
            bottom_margin = 3Plots.mm)
savefig(fig4, joinpath(OUT, "figure_6.pdf"))
println("wrote figures/figure_1.pdf and figures/figure_6.pdf")
