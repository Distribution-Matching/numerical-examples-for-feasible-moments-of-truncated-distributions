#=
figure_2.jl

Draws Figure 2, figures/figure_2.pdf, and prints the parent parameters
quoted in its caption footnote.

For I = [0,1] and mean 0.3, DistributionsFactories.jl finds the location μ
and scale σ at which the truncation to [0,1] of the normal, logistic and
Laplace parents has standard deviation q·sqrt(v_I(0.3)), where v_I(0.3) is
the variance of the truncated exponential with mean 0.3, the largest
possible for a log-concave density (Theorem 1). Panel (a) is q = 0.85 and
panel (b) is q = 0.95, each against the truncated exponential (black).
=#

using DistributionsFactories, Distributions, Plots, LaTeXStrings, Printf
gr()

const MBAR = 0.3

# Truncated exponential ∝ exp(-λx) on [0,1]: mean, and the rate with mean MBAR.
mexp(λ) = 1 / λ - 1 / expm1(λ)
vexp(λ) = 1 / λ^2 - 1 / (4sinh(λ / 2)^2)
function rate_for_mean(m)
    lo, hi = 1e-8, 1e3                          # m < 1/2, so the rate is positive
    for _ in 1:200
        mid = (lo + hi) / 2
        mexp(mid) > m ? (lo = mid) : (hi = mid)
    end
    return (lo + hi) / 2
end

const λ = rate_for_mean(MBAR)
const SDMAX = sqrt(vexp(λ))
@printf("mean %.2f: lambda_I = %.4f, sqrt(v_I) = %.5f\n", MBAR, λ, SDMAX)

const FAMS = ((Normal, "normal", RGB(0.15, 0.25, 0.70)),
              (Logistic, "logistic", RGB(0.90, 0.50, 0.05)),
              (Laplace, "Laplace", RGB(0.10, 0.55, 0.20)))

matched(D, q) = make_dist(truncated(D(); lower = 0, upper = 1), mean = MBAR, std = q * SDMAX)

println("q     normal (mu, sigma)    logistic (mu, sigma)   Laplace (mu, sigma)")
for q in (0.85, 0.95, 0.99)
    @printf("%.2f", q)
    for (D, _, _) in FAMS
        p = matched(D, q).untruncated
        @printf("   (%7.3f, %6.3f)", location(p), scale(p))
    end
    println()
end

function panel(q, title)
    p = plot(; fontfamily = "Computer Modern", framestyle = :box, xlims = (0, 1),
             ylims = (0, 5), xlabel = L"x", ylabel = "density", title = title,
             titlefontsize = 10, guidefontsize = 10, tickfontsize = 8,
             legendfontsize = 7, legend = :topright, foreground_color_legend = :gray70,
             background_color_legend = :white, gridalpha = 0.08)
    for (D, label, col) in FAMS
        d = matched(D, q)
        μ = location(d.untruncated)
        xs = sort(unique(vcat(collect(range(0, 1; length = 401)), 0 < μ < 1 ? [μ] : Float64[])))
        plot!(p, xs, pdf.(d, xs); color = col, lw = 1.4, label = label)
    end
    xs = range(0, 1; length = 401)
    plot!(p, xs, λ .* exp.(-λ .* xs) ./ (-expm1(-λ)); color = :black, lw = 2.2,
          label = "exponential (reference)")
    return p
end

fig = plot(panel(0.85, "(a)"), panel(0.95, "(b)"); layout = (1, 2), size = (620, 205),
           left_margin = 2Plots.mm, bottom_margin = 3Plots.mm)
savefig(fig, joinpath(@__DIR__, "..", "figures", "figure_2.pdf"))
println("wrote figures/figure_2.pdf")
