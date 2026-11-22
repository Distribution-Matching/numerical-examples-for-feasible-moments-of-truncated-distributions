#=
section_5_variance_peaks.jl

Numbers in the remark after Proposition 2 (Section 5): on I = [0,1], along
the ray -μ/σ² = λ with λ < 0 (limiting mean above 1/2), the variance of the
truncation ∝ exp(-λx - x²/(2σ²)) need not increase with σ. For each λ the
script prints the limiting mean, whether the variance is monotone in σ on a
fine grid, the location σ and size (relative to v_I) of its peak, and the
range σ < 1/sqrt(2|λ|) where the proof of Proposition 2 (e)(i) still applies.

Run:  julia --project=. scripts/section_5_variance_peaks.jl
=#

using QuadGK, Printf

function ray_moments(λ, σ)
    lw(x) = -λ * x - x^2 / (2σ^2)
    c = max(lw(0.0), lw(1.0), lw(clamp(-λ * σ^2, 0, 1)))
    g(x) = exp(lw(x) - c)
    Z = quadgk(g, 0, 1; rtol = 1e-12)[1]
    m = quadgk(x -> x * g(x), 0, 1; rtol = 1e-12)[1] / Z
    v = quadgk(x -> (x - m)^2 * g(x), 0, 1; rtol = 1e-12)[1] / Z
    return m, v
end

# truncated exponential with rate λ on [0,1]: mean and variance
texp(λ) = (1 / λ - 1 / expm1(λ), 1 / λ^2 - exp(λ) / expm1(λ)^2)

σs = exp10.(range(-2, 3; length = 20001))
println("λ       limiting mean   monotone   peak σ   peak var / v_I   guaranteed σ <")
for λ in (-0.5, -0.75, -1.0, -2.0, -5.0)
    mlim, vlim = texp(λ)
    vs = [ray_moments(λ, σ)[2] for σ in σs]
    mono = all(diff(vs) .> 0) && vs[end] < vlim
    k = argmax(vs)
    @printf("%-6.2f  %.3f           %-5s      %s   %s           %.2f\n", λ, mlim, mono,
            mono ? "  -  " : @sprintf("%.2f", σs[k]),
            mono ? "  -  " : @sprintf("%.3f", vs[k] / vlim), 1 / sqrt(2abs(λ)))
end
