#=
section_7_ifr_ceiling.jl

Numbers in Section 7 on distributions with increasing failure rate on
[0,1]. At fixed mean the variance of such a distribution is at most that of
the exponential censored at 1, min(E, 1) with E exponential of rate c. The
script prints the peak of this ceiling in mean and standard deviation
coordinates, and a few values against Γ_I and the Bhatia-Davis bound.

Run:  julia --project=. scripts/section_7_ifr_ceiling.jl
=#

using Printf

# censored exponential min(Exp(c), 1): mean, second moment, variance
cm(c) = c < 1e-8 ? 1 - c / 2 : -expm1(-c) / c
cEX2(c) = c < 1e-4 ? 1 - 2c / 3 + c^2 / 4 : 2 * (1 - exp(-c) * (1 + c)) / c^2
cvar(c) = cEX2(c) - cm(c)^2
# truncated exponential on [0,1] with rate λ
tm(λ) = abs(λ) < 1e-6 ? 0.5 - λ / 12 : 1 / λ - 1 / expm1(λ)
tv(λ) = abs(λ) < 1e-4 ? 1 / 12 - λ^2 / 720 : 1 / λ^2 - 1 / (4sinh(λ / 2)^2)

function invert(f, m, lo, hi)          # f decreasing on [lo, hi]
    for _ in 1:200
        mid = (lo + hi) / 2
        f(mid) > m ? (lo = mid) : (hi = mid)
    end
    return (lo + hi) / 2
end
vC(m) = cvar(invert(cm, m, 0.0, 1e4))
vI(m) = tv(invert(tm, m, -1e3, 1e3))

ms = range(1e-4, 1 - 1e-4; length = 20001)
k = argmax(vC.(ms))
@printf("IFR ceiling peaks at sd %.4f, mean %.4f\n", sqrt(vC(ms[k])), ms[k])
println("mean   sd on Γ_I   IFR ceiling   Bhatia-Davis")
for m in (0.1, 0.3, 0.5, 0.7, 0.9)
    @printf("%.1f    %.4f      %.4f        %.4f\n", m, sqrt(vI(m)), sqrt(vC(m)), sqrt(m * (1 - m)))
end
