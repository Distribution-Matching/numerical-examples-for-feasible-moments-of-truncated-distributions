#=
section_5_inner_boundary.jl

Numbers for the inner boundary v_I^in (Section 5, Proposition inner) and for
the open problem on other parents (Section 7).

(i)   v_I^in on I = [0,1] at a few means, against v_I, and its slope at 0.
(ii)  The Swiss data of Section 6: S_n^2 / v_I^in(X̄_n), which is below 1
      exactly when the fitted truncated normal is inner (Corollary mle-inner).
(iii) For I = [0,∞) and several log-concave parents ρ, the coefficient of
      variation r(β) of the truncation of ρ_{μ,σ}, with β = -μ/σ, on a grid,
      a check that it is nondecreasing, and r(0).

Run:  julia --project=. scripts/section_5_inner_boundary.jl
=#

using QuadGK, Roots, Printf

# ---------- (i) v_I^in on [0,1] ----------

# Mean and variance of the truncation to [0,1] of the normal (μ, σ²).
function tn_moments(μ, σ)
    lw(x) = -(x - μ)^2 / (2σ^2)
    c = maximum(lw.((0.0, 1.0, clamp(μ, 0, 1))))
    g(x) = exp(lw(x) - c)
    Z = quadgk(g, 0, 1; rtol = 1e-12)[1]
    m = quadgk(x -> x * g(x), 0, 1; rtol = 1e-12)[1] / Z
    v = quadgk(x -> (x - m)^2 * g(x), 0, 1; rtol = 1e-12)[1] / Z
    return m, v
end

# Truncated exponential on [0,1] with mean m: rate λ and variance v_I(m).
function exp_var(m)
    abs(m - 0.5) < 1e-12 && return 1 / 12
    mean_rate(λ) = abs(λ) < 1e-8 ? 0.5 - λ / 12 : 1 / λ - 1 / expm1(λ)
    λ = find_zero(λ -> mean_rate(λ) - m, (-2000.0, 2000.0), Bisection())
    abs(λ) < 1e-6 && return 1 / 12
    return 1 / λ^2 - exp(λ) / expm1(λ)^2
end

function vin(m)
    m > 0.5 && return vin(1 - m)
    m == 0.5 && return 1 / 12
    σ = find_zero(s -> tn_moments(0.0, s)[1] - m, (1e-3, 1e4), Bisection())
    return tn_moments(0.0, σ)[2]
end

println("(i) inner boundary on [0,1]")
for m in (0.05, 0.1, 0.25, 0.4, 0.5)
    vi, ve = vin(m), exp_var(m)
    @printf("  mean %.2f   sd on Γ^in %.4f   sd on Γ %.4f   CV on Γ^in %.4f\n",
            m, sqrt(vi), sqrt(ve), sqrt(vi) / m)
end
@printf("  half-line CV sqrt(π/2-1) = %.4f\n", sqrt(π / 2 - 1))

# ---------- (ii) Swiss data ----------
# Sample means and biased variances of the four swiss percentages / 100
# (computed in R from datasets::swiss, n = 47).
swiss = [("Agriculture", 0.5065957, 0.0504825),
         ("Examination", 0.1648936, 0.006229244),
         ("Education",   0.1097872, 0.009048891),
         ("Catholic",    0.4114383, 0.1702288)]
println("\n(ii) Swiss data, S^2 / v^in(X̄) and S^2 / v(X̄)")
for (name, m, s2) in swiss
    @printf("  %-12s %.3f   %.3f\n", name, s2 / vin(m), s2 / exp_var(m))
end

# ---------- (iii) other parents on [0,∞) ----------
# Log-densities ℓ with mode at 0; the truncation of ρ_{μ,σ} to [0,∞), scaled by
# 1/σ, has density ∝ exp(ℓ(x + β)) on x ≥ 0 with β = -μ/σ.
log1pexp(z) = z > 0 ? z + log1p(exp(-z)) : log1p(exp(z))
parents = [
    ("normal",                     z -> -z^2 / 2),
    ("logistic",                   z -> -z - 2log1pexp(-z)),
    ("Gumbel",                     z -> -(z + exp(-z))),
    ("reflected Gumbel",           z -> z - exp(z)),
    ("generalised normal α=1.5",   z -> -abs(z)^1.5),
    ("generalised normal α=3",     z -> -abs(z)^3),
    ("generalised normal α=4",     z -> -z^4),
    ("Laplace",                    z -> -abs(z)),
]

function cv(ℓ, β)
    c = ℓ(max(β, 0.0))                    # log-density at the mode on [0,∞)
    g(x) = exp(ℓ(x + β) - c)
    # scale of the bulk: up to where the log-density has dropped by 60
    U = find_zero(x -> ℓ(x + β) - c + 60, (max(-β, 0.0), 1e6), Bisection())
    Z = quadgk(g, 0, U / 8, U; rtol = 1e-11)[1]
    m = quadgk(x -> x * g(x), 0, U / 8, U; rtol = 1e-11)[1] / Z
    v = quadgk(x -> (x - m)^2 * g(x), 0, U / 8, U; rtol = 1e-11)[1] / Z
    return sqrt(v) / m
end

println("\n(iii) CV r(β) on [0,∞), β ∈ [-4, 6]")
βs = range(-4, 6; length = 401)
for (name, ℓ) in parents
    r = [cv(ℓ, β) for β in βs]
    mono = all(diff(r) .> -1e-9)
    @printf("  %-26s r(0) = %.4f   nondecreasing: %s   r(6) = %.4f\n",
            name, cv(ℓ, 0.0), mono, r[end])
end
