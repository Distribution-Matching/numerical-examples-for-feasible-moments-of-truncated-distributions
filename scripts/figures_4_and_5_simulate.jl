#=
figures_4_and_5_simulate.jl

Generates the data for Figures 4 and 5 (about 3 minutes with 8 threads;
run as `julia -t auto --project=. scripts/figures_4_and_5_simulate.jl`).

For I = [0,1], every point (m, σ̄) with 0 < σ̄² < m(1-m) is the mean and
standard deviation of a unique density f_θ(x) ∝ exp(θ₁x + θ₂x²) on [0,1]
(Brown 1986, Theorem 3.6): a truncated normal (θ₂ < 0) below Γ_I, a
truncated exponential (θ₂ = 0) on Γ_I, and a density with convex logarithm
(θ₂ > 0, an "upside-down Gaussian") above Γ_I. For each grid point we draw
R = 10,000 samples of size n from f_θ, with common random numbers (the same
R × n uniforms at every point, transformed by the quantile function), and
record the fraction whose sample mean and
sample variance (divisor n) lie strictly below Γ_I, i.e. the probability
that the maximum likelihood fit is a truncated normal.

The picture is symmetric under x ↦ 1-x, so only m ≤ 1/2 is simulated and
the result is mirrored.

For I = [0,∞) the model has no member above Γ_I (the line σ̄ = m), so
only truncated normal data are simulated, and we record the probability that
s ≥ m, that is, that the maximum likelihood estimator does not exist (the
fit is the exponential on Γ_I). By scale invariance this depends only on the
coefficient of variation c = σ̄/m < 1, which is gridded.

Writes data/mle_regime_unit_n{N}.dat (columns m, sd, p; one
scanline per m, blank-line separated, NaN outside the region),
data/mle_regime_half.dat (columns c and P(s ≥ m) for each n) and
data/mle_regime_curves.dat (m, sqrt(v_I(m)), sqrt(m(1-m))).
=#

using Printf, Random, Distributions, DistributionsFactories

const K = 4000                                   # Simpson panels (even)
const XS = collect(range(0.0, 1.0; length = K + 1))
const WS = [i == 0 || i == K ? 1.0 : (isodd(i) ? 4.0 : 2.0) for i in 0:K] ./ (3K)

# Normalised weights of f_θ on the grid, computed stably.
function weights(θ)
    lp = θ[1] .* XS .+ θ[2] .* XS .^ 2
    c = maximum(lp)
    w = exp.(lp .- c) .* WS
    Z = sum(w)
    return w ./ Z, log(Z) + c
end

# Solve ∇ψ(θ) = t by damped Newton on the convex function ψ(θ) - θ·t.
function solve_theta(m, v)
    t = (m, v + m^2)
    θ = [m / v, -1 / (2v)]                       # Gaussian guess
    obj(θ) = weights(θ)[2] - θ[1] * t[1] - θ[2] * t[2]
    if !isfinite(obj(θ)) || v > 0.9 * m * (1 - m)
        θ = [0.0, 0.0]
    end
    for _ in 1:200
        w, _ = weights(θ)
        e1 = sum(w .* XS); e2 = sum(w .* XS .^ 2)
        e3 = sum(w .* XS .^ 3); e4 = sum(w .* XS .^ 4)
        g = [e1 - t[1], e2 - t[2]]
        H = [e2-e1^2 e3-e1*e2; e3-e1*e2 e4-e2^2]
        maximum(abs.(g)) < 1e-9 && return θ
        d = -(H \ g)
        f0 = obj(θ); s = 1.0
        while obj(θ .+ s .* d) > f0 + 1e-4 * s * (g' * d) && s > 1e-12
            s /= 2
        end
        θ = θ .+ s .* d
    end
    # The line search can stall at round-off level; accept if close.
    w, _ = weights(θ)
    r = max(abs(sum(w .* XS) - t[1]), abs(sum(w .* XS .^ 2) - t[2]))
    r < 1e-7 || @warn "Newton did not converge" m v r
    return θ
end

# Boundary v_I(m), tabulated through the rate λ and interpolated in m.
mexp(λ) = abs(λ) < 1e-4 ? 0.5 - λ / 12 : 1 / λ - 1 / expm1(λ)
vexp(λ) = abs(λ) < 1e-3 ? 1 / 12 - λ^2 / 720 : 1 / λ^2 - 1 / (4sinh(λ / 2)^2)
const LAM = collect(range(3000.0, -3000.0; length = 600_001))   # m increasing
const MTAB = mexp.(LAM)
const VTAB = vexp.(LAM)
function vI(m)
    m <= MTAB[1] && return m^2
    m >= MTAB[end] && return (1 - m)^2
    j = searchsortedlast(MTAB, m)
    a = (m - MTAB[j]) / (MTAB[j+1] - MTAB[j])
    return (1 - a) * VTAB[j] + a * VTAB[j+1]
end

# Common random numbers: for each n, one fixed R × n matrix of uniforms is
# used at every grid point, each density being sampled by inverse transform.
const NS = (5, 10, 20, 100)
const R = 10_000
const UNIF = Dict(n => rand(MersenneTwister(20261110 + n), R, n) for n in NS)

# Fraction of the R samples of size n from f_θ with sample point below Γ_I.
function p_below(θ, n)
    U = UNIF[n]
    w, _ = weights(θ)
    F = cumsum(w); F ./= F[end]
    below = 0
    for r in 1:R
        s1 = 0.0; s2 = 0.0
        for k in 1:n
            u = U[r, k]
            j = min(searchsortedfirst(F, u), length(F))
            x = j == 1 ? 0.0 :
                XS[j-1] + (u - F[j-1]) / (F[j] - F[j-1]) * (XS[j] - XS[j-1])
            s1 += x; s2 += x^2
        end
        mm = s1 / n; ss = s2 / n - mm^2
        below += ss < vI(mm)
    end
    return below / R
end

const NM = 200                                   # cells in m (even)
const NSD = 160                                  # cells in sd
const MS = [(i - 0.5) / NM for i in 1:NM]
const SDS = [(j - 0.5) * 0.5 / NSD for j in 1:NSD]

function main()
    P = Dict(n => fill(NaN, NM, NSD) for n in NS)
    cells = [(i, j) for i in 1:NM÷2, j in 1:NSD if SDS[j]^2 < 0.97 * MS[i] * (1 - MS[i])]
    Threads.@threads for (i, j) in cells
        θ = solve_theta(MS[i], SDS[j]^2)
        for n in NS
            p = p_below(θ, n)
            P[n][i, j] = p
            P[n][NM+1-i, j] = p                  # mirror x ↦ 1-x
        end
    end
    dir = joinpath(@__DIR__, "..", "data")
    for n in NS
        open(joinpath(dir, "mle_regime_unit_n$(n).dat"), "w") do io
            println(io, "m sd p")
            for i in 1:NM
                for j in 1:NSD
                    @printf(io, "%.5f %.5f %s\n", MS[i], SDS[j],
                            isnan(P[n][i, j]) ? "nan" : @sprintf("%.4f", P[n][i, j]))
                end
                println(io)
            end
        end
    end
    half(dir)
    open(joinpath(dir, "mle_regime_curves.dat"), "w") do io
        println(io, "m gamma outer")
        for m in range(0.0, 1.0; length = 401)
            @printf(io, "%.5f %.6f %.6f\n", m, sqrt(vI(m)), sqrt(m * (1 - m)))
        end
    end
end

# I = [0,∞): the quantile function of the truncated normal on [0,∞) with
# mean 1 and coefficient of variation c < 1, constructed with
# DistributionsFactories.jl. It is inverted through the log upper tail, so
# that the very negative parent locations needed as c → 1 stay accurate.
function half_quantile(c)
    p = make_dist(truncated(Normal(); lower = 0), mean = 1.0, std = c).untruncated
    μ, σ = p.μ, p.σ
    L = logccdf(Normal(), -μ / σ)
    return u -> μ + σ * invlogccdf(Normal(), L + log1p(-u))
end

# Fraction of the R samples of size n with s ≥ m.
function p_fail(c, n)
    U = UNIF[n]
    q = half_quantile(c)
    fail = 0
    x = zeros(n)
    for r in 1:R
        for k in 1:n
            x[k] = q(U[r, k])
        end
        mm = sum(x) / n
        fail += sum(abs2, x .- mm) / n >= mm^2
    end
    return fail / R
end

const CS = collect(range(0, 1; length = 202)[2:end-1])   # 200 values of c in (0,1)

function half(dir)
    P = zeros(length(CS), length(NS))
    Threads.@threads for k in eachindex(CS)
        for (j, n) in enumerate(NS)
            P[k, j] = p_fail(CS[k], n)
        end
    end
    open(joinpath(dir, "mle_regime_half.dat"), "w") do io
        println(io, "c " * join(("n$n" for n in NS), " "))
        for k in eachindex(CS)
            @printf(io, "%.4f %s\n", CS[k], join((@sprintf("%.4f", P[k, j]) for j in eachindex(NS)), " "))
        end
    end
end

main()
