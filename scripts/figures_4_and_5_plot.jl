#=
figures_4_and_5_plot.jl

Draws Figures 4 and 5 from the data written by figures_4_and_5_simulate.jl.
figures/figure_4.pdf (Figure 4), I = [0,1], n = 5, 10, 20, 100: the colour
is the probability that the sample point falls on the other side of Γ_I from
the true point, so that the maximum likelihood fit is of the other type:
P(below) for true points above Γ_I (purple) and P(above) for true points
below it (red), with the two colour scales side by side on the right. The
green curve is the level 1/2 of P(below).
figures/figure_5.pdf (Figure 5), I = [0,∞): for truncated normal data, the
probability that S_n ≥ X̄_n, so that the maximum likelihood estimator does not
exist, against the coefficient of variation σ̄/μ̄, for the same n.
=#

using Plots, LaTeXStrings, Printf
gr()

const DIR = joinpath(@__DIR__, "..", "data")
const NS = (5, 10, 20, 100)
# Colour stops, placed towards 0 so that probabilities of 0.1-0.2 are
# already clearly coloured.
const STOPS = [0.0, 0.04, 0.15, 0.4, 1.0]
const PCOLS = [RGB(1, 1, 1), RGB(0.86, 0.81, 0.93), RGB(0.68, 0.57, 0.82),
               RGB(0.42, 0.26, 0.60), RGB(0.20, 0.06, 0.35)]   # true above, sample below
const RCOLS = [RGB(1, 1, 1), RGB(0.99, 0.85, 0.80), RGB(0.96, 0.60, 0.52),
               RGB(0.82, 0.25, 0.20), RGB(0.50, 0.03, 0.08)]   # true below, sample above
const PURPLE = cgrad(PCOLS, STOPS)
const RED = cgrad(RCOLS, STOPS)
const BLUE = RGB(0.16, 0.24, 0.55)             # Γ_I, as in Figure 1
const GREEN = RGB(0.05, 0.70, 0.25)            # level 1/2 of P(below)
const EMAX = 1.0                                # top of the colour scale

function readgrid(file)
    rows = filter(!isempty, [split(l) for l in eachline(file)][2:end])
    ms = sort(unique(parse.(Float64, first.(rows))))
    sds = sort(unique(parse.(Float64, getindex.(rows, 2))))
    P = fill(NaN, length(sds), length(ms))
    for r in rows
        i = searchsortedfirst(ms, parse(Float64, r[1]))
        j = searchsortedfirst(sds, parse(Float64, r[2]))
        P[j, i] = r[3] == "nan" ? NaN : parse(Float64, r[3])
    end
    return ms, sds, P
end

curves = [split(l) for l in eachline(joinpath(DIR, "mle_regime_curves.dat"))][2:end]
cm = parse.(Float64, first.(curves))
cgam = parse.(Float64, getindex.(curves, 2))
cout = parse.(Float64, getindex.(curves, 3))

gammaI(m) = (k = clamp(searchsortedlast(cm, m), 1, length(cm) - 1);
             t = (m - cm[k]) / (cm[k+1] - cm[k]); (1 - t) * cgam[k] + t * cgam[k+1])

# Probability of the wrong regime, given P(below), on one true side only
# (NaN elsewhere): side = :above gives P(below), side = :below gives P(above).
wrong(P, ms, sds, ceiling, side) =
    [isnan(P[j, i]) || ((sds[j] < ceiling(ms[i])) != (side == :below)) ? NaN :
     (side == :below ? 1 - P[j, i] : P[j, i]) for j in eachindex(sds), i in eachindex(ms)]
# Both sides in one heatmap: -P(below) above Γ_I, +P(above) below it, drawn
# with a colour map running purple (-EMAX), white (0), red (+EMAX).
const SIGNED = cgrad(vcat(reverse(PCOLS), RCOLS[2:end]),
                     vcat(0.5 .- reverse(STOPS) ./ 2, 0.5 .+ STOPS[2:end] ./ 2))
function twoheat(ms, sds, P, ceiling; kw...)
    A = wrong(P, ms, sds, ceiling, :above)
    B = wrong(P, ms, sds, ceiling, :below)
    W = map((a, b) -> isnan(a) ? b : -a, A, B)
    return heatmap(ms, sds, W; c = SIGNED, clims = (-EMAX, EMAX), kw...)
end
# Level 1/2 of P(below). Cells outside the region (NaN) take the last valid
# value below them in their column, so the contour is not drawn along the
# outer boundary.
function half_contour!(p, ms, sds, P)
    Q = copy(P)
    for i in axes(Q, 2), j in axes(Q, 1)
        isnan(Q[j, i]) && (Q[j, i] = j > 1 ? Q[j-1, i] : 1.0)
    end
    contour!(p, ms, sds, Q; levels = [0.5], color = GREEN, lw = 1.5,
             colorbar = false)
end

common = (fontfamily = "Computer Modern", colorbar = false,
          framestyle = :box, xlims = (0, 1), legend = false, grid = false,
          xlabel = L"\bar\mu", ylabel = L"\bar\sigma", titlefontsize = 10,
          guidefontsize = 10, tickfontsize = 8, xticks = 0:0.2:1)

frame!(p, xmax, ymax) = plot!(p, [0, xmax, xmax, 0, 0], [0, 0, ymax, ymax, 0];
                              color = :black, lw = 0.8)
const HMAX = 5.0          # half-line axes; the picture is scale invariant

panels = []
for n in NS
    ms, sds, P = readgrid(joinpath(DIR, "mle_regime_unit_n$(n).dat"))
    p = twoheat(ms, sds, P, gammaI; common..., ylims = (0, 0.5), yticks = 0:0.1:0.5,
                title = L"n=%$n")
    plot!(p, cm, cgam; color = BLUE, lw = 2.0)
    half_contour!(p, ms, sds, P)
    # White mask outside the region, hiding the contour's continuation there.
    plot!(p, vcat(cm, 1, 0), vcat(cout, 0.5, 0.5); seriestype = :shape,
          fillcolor = :white, linecolor = :white, lw = 0)
    plot!(p, cm, cout; color = :gray45, lw = 0.9)
    frame!(p, 1.0, 0.5)
    push!(panels, p)
end

# Two parallel vertical colour scales, purple then red, ticks on the right.
vbar(cmap, ticks; right = true) = heatmap([0, 1], range(0, EMAX; length = 256),
                            repeat(collect(range(0, EMAX; length = 256)), 1, 2);
                            c = cmap, clims = (0, EMAX), colorbar = false,
                            xticks = false, yticks = ticks ? (0:0.25:1) : false,
                            xlims = (-0.5, 1.5), ylims = (0, EMAX),
                            ymirror = right, framestyle = :box,
                            fontfamily = "Computer Modern", tickfontsize = 8,
                            left_margin = 0Plots.mm, right_margin = 0Plots.mm)
blank() = plot(; framestyle = :none, grid = false, ticks = false, legend = false)
# The green level-1/2 curve of the panels is marked on the purple scale.
pbar = vbar(PURPLE, false; right = false)
plot!(pbar; yticks = ([0, 0.5, 1], ["0", "0.5", "1"]))   # the 1/2 level is labelled
plot!(pbar, [-0.5, 1.5], [0.5, 0.5]; color = GREEN, lw = 3.5, legend = false)
scales = plot(blank(), pbar, vbar(RED, true), blank();
              layout = @layout([a{0.12h}; [b c]; d{0.12h}]))

fig = plot(plot(panels...; layout = (2, 2)), scales;
           layout = @layout([a b{0.13w}]), size = (640, 400),
           left_margin = 2Plots.mm, bottom_margin = 2Plots.mm, right_margin = 3Plots.mm)
savefig(fig, joinpath(@__DIR__, "..", "figures", "figure_4.pdf"))

# Figure 5: half-line, probability that the estimator does not exist, in
# reds darkening with n.
const HCOLS = (RGB(0.97, 0.58, 0.48), RGB(0.88, 0.38, 0.30), RGB(0.74, 0.16, 0.14),
               RGB(0.45, 0.02, 0.07))
rows = [split(l) for l in eachline(joinpath(DIR, "mle_regime_half.dat"))][2:end]
cs = parse.(Float64, first.(rows))
hp = plot(; xlims = (0.5, 1), ylims = (0, 0.45), xticks = 0.5:0.1:1, framestyle = :box, gridalpha = 0.08,
          xlabel = L"\bar\sigma/\bar\mu", ylabel = "Probability that " * L"S_n\geq\bar X_n",
          legend = :topleft, fontfamily = "Computer Modern", guidefontsize = 10,
          tickfontsize = 8, legendfontsize = 8, background_color_legend = :white,
          foreground_color_legend = :gray70, size = (460, 230),
          left_margin = 3Plots.mm, bottom_margin = 3Plots.mm)
for (j, n) in enumerate(NS)
    plot!(hp, cs, parse.(Float64, getindex.(rows, j + 1)); color = HCOLS[j], lw = 1.8,
          label = L"n=%$n")
end
savefig(hp, joinpath(@__DIR__, "..", "figures", "figure_5.pdf"))
println("wrote figures/figure_4.pdf and figures/figure_5.pdf")
