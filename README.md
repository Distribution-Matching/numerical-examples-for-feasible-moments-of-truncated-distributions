# Numerical examples for feasible moments of truncated distributions

Julia scripts that reproduce the figures of the paper
*Feasible moments of singly or doubly truncated Gaussians and other
log-concave distributions* (Moka & Nazarathy).

The generated figures are included in [`figures/`](figures), and the
simulation output behind Figures 4 and 5 in [`data/`](data).

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
scripts/run_all.jl                   # all of the above
figures/figure_N.pdf                 # the figures
data/                                # simulation output for Figures 4 and 5
```

## Running

With Julia 1.10 or later, from the repository root:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia -t auto --project=. scripts/run_all.jl
```

The simulation for Figures 4 and 5 takes about 3 minutes with 8 threads;
the other figures take seconds. Each script can also be run on its own,
for example `julia --project=. scripts/figure_2.jl`. The plotting script
for Figures 4 and 5 reads the data written by the simulation script.
