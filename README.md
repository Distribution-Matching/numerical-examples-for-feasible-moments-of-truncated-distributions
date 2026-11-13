# Numerical examples for feasible moments of truncated distributions

Julia scripts that reproduce the figures and Table 1 of the paper
*Feasible means and variances of truncated log-concave distributions* (Moka & Nazarathy).

The generated figures are included in [`figures/`](figures), and the
simulation output behind Figures 4 and 5 and Table 1 in [`data/`](data).

## Companion artefacts

- Julia package used here for moment matching:
  <https://github.com/Distribution-Matching/DistributionsFactories.jl>
- R and Python counterparts:
  <https://github.com/Distribution-Matching/distsfactory-r>,
  <https://github.com/Distribution-Matching/distsfactory-python>

## Layout

```
scripts/figures_1_and_6.jl           # Figure 1: feasible region for I = [0,1] and [0,∞)
                                     # Figure 6: conjectured Student-t boundaries
scripts/figure_2.jl                  # Figure 2: normal, logistic and Laplace truncations
                                     #   matched to a mean and standard deviation with
                                     #   DistributionsFactories.make_dist; also prints
                                     #   the parent parameters quoted in the paper
scripts/figure_3.jl                  # Figure 3: Gaussian rays μ = -λσ² and their approach
                                     #   to the boundary (Proposition 2)
scripts/figures_4_and_5_simulate.jl  # Figures 4 and 5: Monte Carlo study of the maximum
                                     #   likelihood regime (writes data/)
scripts/figures_4_and_5_plot.jl      # Figures 4 and 5: plots from data/
scripts/table_1.jl                   # Table 1: size and power of the diagnostic test of
                                     #   Section 6 for I = [0,1], with 100,000 replications
                                     #   per cell, and the first-order (CLT) predictions
                                     #   (writes data/table_1.dat)
scripts/table_1_size_check.jl        # size of the half-line test at the exponential for
                                     #   n up to 5000, quoted in the discussion of Table 1
scripts/run_all.jl                   # all of the above except the size check
figures/figure_N.pdf                 # the figures
data/                                # simulation output for Figures 4 and 5 and Table 1
```

## Running

With Julia 1.10 or later, from the repository root:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia -t auto --project=. scripts/run_all.jl
```

The simulation for Figures 4 and 5 takes about 3 minutes with 8 threads,
Table 1 about a minute, and the other figures seconds. Each script can also be run on its own,
for example `julia --project=. scripts/figure_2.jl`. The plotting script
for Figures 4 and 5 reads the data written by the simulation script.
