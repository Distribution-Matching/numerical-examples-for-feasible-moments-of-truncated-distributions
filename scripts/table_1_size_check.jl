#=
table_1_size_check.jl

Size of the diagnostic test of Section 6 at the exponential on I = [0,∞),
for larger n than in Table 1 (200,000 replications per n, two seeds). On
the half-line v_I(m) = m^2 and ω(m)^2 = 4 m^4, so the test rejects when
sqrt(n) (S_n^2 - X̄_n^2) > Φ⁻¹(0.95) 2 X̄_n^2. Prints the sizes, which stay
near 0.053 for n between 500 and 5000 (skewness of the summands, see the
discussion of Table 1).
=#

using Random, Distributions, Printf

const c = quantile(Normal(), 0.95)

function size_exp(n, R, seed)
    rng = MersenneTwister(seed); x = zeros(n); rej = 0
    for _ in 1:R
        for k in 1:n
            x[k] = randexp(rng)
        end
        m = sum(x) / n
        s2 = sum(abs2, x .- m) / n
        rej += sqrt(n) * (s2 - m^2) > c * 2m^2
    end
    return rej / R
end

for n in (100, 200, 500, 1000, 2000, 5000)
    @printf("n = %5d   size = %.4f  %.4f\n", n, size_exp(n, 200_000, 7 + n), size_exp(n, 200_000, 99 + n))
end
