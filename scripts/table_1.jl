#=
table_1.jl

Finite-sample calibration of the diagnostic test of Section 6 for I = [0,1].
With T_n = sqrt(n) G_n / ω(X̄_n), three versions of the test at α = 0.05:

  asymptotic       reject when T_n > Φ⁻¹(1-α), the test (14);
  bias-corrected   G_n is replaced by G_n + S_n² ω(X̄_n)² / (2 n v_I(X̄_n)²),
                   which removes the leading bias -σ̄²(1 + v_I''/2)/n, since
                   1 + v_I''/2 = ω²/(2 v_I²) on Γ_I;
  tabulated        reject when T_n > q_n(X̄_n), where q_n(μ) is the
                   (1-α)-quantile of T_n under the truncated exponential
                   ξ_{I,μ}, tabulated on a grid of μ in (0, 1/2] and used at
                   min(X̄_n, 1 - X̄_n) by the reflection symmetry. This is the
                   parametric bootstrap at ξ_{I,X̄_n}, with the bootstrap
                   quantiles computed once in advance.

The quantile table uses 40,000 replications per grid point, and each
rejection rate 100,000 replications (independent seeds). Prints the rows
and writes data/table_1.dat (takes several minutes). The asymptotic and
tabulated rows give Table 1 of the paper.
=#

include(joinpath(@__DIR__, "table_1_core.jl"))

const GRID = vcat(0.005, 0.01:0.01:0.5)
const RQ = 40_000
const RC = 100_000
const α = 0.05
const zα = quantile(Normal(), 1 - α)

stat(m, s2, n) = sqrt(n) * (s2 - vI01(m)) / sqrt(ω201(m))

function sample_moments!(x, rng, d)
    for k in eachindex(x)
        x[k] = rand(rng, d)
    end
    n = length(x); m = sum(x) / n
    return m, sum(abs2, x .- m) / n
end

# q_n(μ) on GRID, from RQ replications under ξ_{I,μ}.
function quantile_table(n; seed = 1)
    q = zeros(length(GRID))
    for (j, μ) in enumerate(GRID)
        rng = MersenneTwister(seed + 1000n + j)
        d = TExp01(rate(μ)); x = zeros(n); T = zeros(RQ)
        for r in 1:RQ
            m, s2 = sample_moments!(x, rng, d)
            T[r] = stat(m, s2, n)
        end
        q[j] = quantile(T, 1 - α)
    end
    return q
end

function qinterp(q, m)
    μ = clamp(min(m, 1 - m), GRID[1], GRID[end])
    j = min(searchsortedlast(GRID, μ), length(GRID) - 1)
    a = (μ - GRID[j]) / (GRID[j+1] - GRID[j])
    return (1 - a) * q[j] + a * q[j+1]
end

function rates(d, n, q; seed = 777)
    rng = MersenneTwister(seed + n); x = zeros(n)
    ra = rb = rt = 0
    for _ in 1:RC
        m, s2 = sample_moments!(x, rng, d)
        T = stat(m, s2, n)
        Tb = T + sqrt(n) * s2 * ω201(m) / (2n * vI01(m)^2) / sqrt(ω201(m))
        ra += T > zα
        rb += Tb > zα
        rt += T > qinterp(q, m)
    end
    return ra / RC, rb / RC, rt / RC
end

function calibration_main()
    dens = [("truncated exponential, mean 0.05", TExp01(rate(0.05))),
            ("truncated exponential, mean 0.25", TExp01(rate(0.25))),
            ("truncated exponential, mean 0.5", TExp01(rate(0.5))),
            ("truncated normal N(-2,1)", truncated(Normal(-2.0, 1.0), 0.0, 1.0)),
            ("Beta(2,2)", Beta(2, 2)),
            ("Beta(0.5,0.5)", Beta(0.5, 0.5)),
            ("Beta(0.8,0.8)", Beta(0.8, 0.8)),
            ("Beta(0.5,1.5)", Beta(0.5, 1.5)),
            ("mixture B(3,12),B(12,3)", MixtureModel([Beta(3, 12), Beta(12, 3)])),
            ("mixture B(5,15),B(15,5)", MixtureModel([Beta(5, 15), Beta(15, 5)]))]
    tabs = Dict(n => quantile_table(n) for n in NS)
    path = joinpath(@__DIR__, "..", "data", "table_1.dat")
    open(path, "w") do io
        println(io, "# I = [0,1], alpha = 0.05: density | method | rejection rate for n = ",
                join(NS, ", "), "   (table RQ = $RQ per grid point, RC = $RC)")
        for (name, d) in dens
            res = [rates(d, n, tabs[n]) for n in NS]
            for (k, meth) in enumerate(("asymptotic", "bias-corrected", "tabulated"))
                line = @sprintf("%-34s %-15s %s", name, meth,
                                join((@sprintf("%.3f", r[k]) for r in res), "  "))
                println(io, line); println(line)
            end
        end
        println(io, "# q_n(mu) at mu = 0.05, 0.25, 0.5: ",
                join((@sprintf("n=%d: %.3f %.3f %.3f", n, qinterp(tabs[n], 0.05),
                               qinterp(tabs[n], 0.25), qinterp(tabs[n], 0.5)) for n in NS), "; "))
    end
end

calibration_main()
