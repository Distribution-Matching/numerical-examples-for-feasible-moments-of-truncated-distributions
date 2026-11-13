#=
table_1.jl

Table 1: simulation study of the diagnostic test of Section 6 (Corollary
cor:boundary-test): reject H_0 (f log-concave on I) when

    S_n^2 / v_I(X̄_n) > 1 + Φ⁻¹(1-α) sqrt((2 + v_I''(X̄_n)) / n),

equivalently sqrt(n) G_n > Φ⁻¹(1-α) ω(X̄_n) with G_n = S_n^2 - v_I(X̄_n) and
ω(m)² = v_I(m)²(2 + v_I''(m)), computed for I = [0,1] through the Langevin
function as ω² = (L''' + 2L'² - L''²/L')/16 at z = L⁻¹(1-2m) (Remark 4).
Here X̄_n and S_n^2 are the sample mean and variance with divisor n.

Each cell uses R = 100,000 replications at α = 0.05, for I = [0,1] (the
half-line case is the classical coefficient-of-variation test; its size at
large n is checked in table_1_size_check.jl). For each density the script also reports the ratio
σ̄²/v_I(μ̄), which locates the true point relative to Γ_I, and the
first-order CLT prediction of the rejection probability,

    1 - Φ( (Φ⁻¹(1-α) ω_ξ - sqrt(n) (σ̄² - v_I(μ̄))) / ω_f ),

with ω_f² = Var((X-μ̄)² - v_I'(μ̄) X) and ω_ξ = ω(μ̄) (Proposition 3).
Prints the rows and writes data/table_1.dat. Takes about a minute.
=#

using Printf, Random, Distributions, QuadGK

setprecision(BigFloat, 256)

# ---------------------------------------------------------------- I = [0,1]
# Mean, variance, ω² of the truncated exponential ∝ e^{-λx} on [0,1],
# through the Langevin function at z = λ/2 (BigFloat avoids cancellation).
function langevin_moments(λ)
    z = BigFloat(λ) / 2
    c = coth(z); s2 = csch(z)^2
    L = c - 1 / z
    L1 = 1 / z^2 - s2
    L2 = -2 / z^3 + 2 * s2 * c
    L3 = 6 / z^4 - 4 * s2 * c^2 - 2 * s2^2
    return Float64((1 - L) / 2), Float64(L1 / 4), Float64((L3 + 2L1^2 - L2^2 / L1) / 16),
           Float64(-L2 / 8 / (L1 / 4))           # v_I'(μ̄) = κ₃/κ₂
end

# Tables in m (increasing as λ decreases); λ = 0 is skipped (removable singularity).
const LAM = [λ for λ in range(2000.0, -2000.0; length = 200_001) if abs(λ) > 1e-7]
const TAB = langevin_moments.(LAM)
const MT = [t[1] for t in TAB]
const VT = [t[2] for t in TAB]
const WT = [t[3] for t in TAB]
const DT = [t[4] for t in TAB]
function interp(m, Y, lo, hi)
    m <= MT[1] && return lo(m)
    m >= MT[end] && return hi(m)
    j = searchsortedlast(MT, m)
    a = (m - MT[j]) / (MT[j+1] - MT[j])
    return (1 - a) * Y[j] + a * Y[j+1]
end
vI01(m) = interp(m, VT, m -> m^2, m -> (1 - m)^2)
ω201(m) = interp(m, WT, m -> 4m^4, m -> 4(1 - m)^4)
dv01(m) = interp(m, DT, m -> 2m, m -> -2(1 - m))

# Rate of the truncated exponential on [0,1] with mean μ, by bisection.
function rate(μ)
    lo, hi = -2000.0, 2000.37                    # hi avoids the midpoint λ = 0
    for _ in 1:200
        mid = (lo + hi) / 2
        langevin_moments(mid)[1] > μ ? (lo = mid) : (hi = mid)
    end
    return (lo + hi) / 2
end

# Truncated exponential on [0,1] with rate λ, sampled by inversion.
struct TExp01 <: ContinuousUnivariateDistribution
    λ::Float64
end
Base.rand(rng::AbstractRNG, d::TExp01) = (u = rand(rng); -log1p(-u * (-expm1(-d.λ))) / d.λ)

# ------------------------------------------------------------- the test
function rejection_rate(d, n, R, vI, ω2; α = 0.05, seed = 20261113)
    rng = MersenneTwister(seed + n)
    c = quantile(Normal(), 1 - α)
    x = zeros(n); rej = 0
    for _ in 1:R
        for k in 1:n
            x[k] = rand(rng, d)
        end
        m = sum(x) / n
        s2 = sum(abs2, x .- m) / n
        rej += sqrt(n) * (s2 - vI(m)) > c * sqrt(ω2(m))
    end
    return rej / R
end

# Moments of f by quadrature (or closed form), for the ratio and prediction.
function moments(d, I01::Bool)
    if d isa TExp01
        μ, v, _, _ = langevin_moments(d.λ)
        k3 = quadgk(x -> (x - μ)^3 * texpdf(d, x), 0, 1)[1]
        k4 = quadgk(x -> (x - μ)^4 * texpdf(d, x), 0, 1)[1]
        return μ, v, k3, k4
    end
    μ = mean(d); v = var(d)
    hi = I01 ? 1.0 : Inf
    m3 = quadgk(x -> (x - μ)^3 * pdf(d, x), 0, hi)[1]
    m4 = quadgk(x -> (x - μ)^4 * pdf(d, x), 0, hi)[1]
    return μ, v, m3, m4
end
texpdf(d::TExp01, x) = d.λ * exp(-d.λ * x) / (-expm1(-d.λ))

function prediction(d, n, vI, dv, ω2, I01; α = 0.05)
    μ, v, m3, m4 = moments(d, I01)
    b = dv(μ)
    # ω_f² = Var((X-μ)² - b X) = m4 - v² - 2 b m3 + b² v
    ωf = sqrt(m4 - v^2 - 2b * m3 + b^2 * v)
    c = quantile(Normal(), 1 - α)
    return v / vI(μ), 1 - cdf(Normal(), (c * sqrt(ω2(μ)) - sqrt(n) * (v - vI(μ))) / ωf)
end

const NS = (20, 50, 100, 500)
const R = 100_000

function main()
    rows = Tuple{String,String,Float64,Vector{Float64},Vector{Float64}}[]
    A = [("truncated exponential, mean 0.05", TExp01(rate(0.05))),
         ("truncated exponential, mean 0.25", TExp01(rate(0.25))),
         ("truncated exponential, mean 0.5", TExp01(rate(0.5))),
         ("truncated normal N(-2,1)", truncated(Normal(-2.0, 1.0), 0.0, 1.0)),
         ("Beta(2,2)", Beta(2, 2)),
         ("Beta(0.5,0.5)", Beta(0.5, 0.5)),
         ("Beta(0.8,0.8)", Beta(0.8, 0.8)),
         ("Beta(0.5,1.5)", Beta(0.5, 1.5)),
         ("mixture B(3,12),B(12,3)", MixtureModel([Beta(3, 12), Beta(12, 3)])),
         ("mixture B(5,15),B(15,5)", MixtureModel([Beta(5, 15), Beta(15, 5)]))]
    for (name, d) in A
        ratio, _ = prediction(d, 20, vI01, dv01, ω201, true)
        p = [rejection_rate(d, n, R, vI01, ω201) for n in NS]
        q = [prediction(d, n, vI01, dv01, ω201, true)[2] for n in NS]
        push!(rows, ("A", name, ratio, p, q))
    end
    path = joinpath(@__DIR__, "..", "data", "table_1.dat")
    open(path, "w") do io
        println(io, "# I = [0,1]: density | sigma^2/v_I(mu) | rejection rate for n = ", join(NS, ", "),
                " | CLT prediction for the same n   (R = $R, alpha = 0.05)")
        for (pan, name, r, p, q) in rows
            line = @sprintf("%-34s %.3f  %s  |  %s", name, r,
                            join((@sprintf("%.3f", x) for x in p), "  "),
                            join((@sprintf("%.3f", x) for x in q), "  "))
            println(io, line); println(line)
        end
    end
end

main()
