# Numerical examples for feasible moments of truncated distributions

Julia and R scripts that reproduce the figures, the table and the
numerical results of the paper
*Feasible means and variances of truncated log-concave distributions* (Moka & Nazarathy).

The generated figures are included in [`figures/`](figures), and the
simulation output behind Figures 4 and 5, Table 1 and Section 6 in [`data/`](data).

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
scripts/table_1.jl                   # Table 1: rejection rates of the diagnostic test of
                                     #   Section 6 for I = [0,1], asymptotic and with
                                     #   tabulated null quantiles (also the bias-corrected
                                     #   version mentioned in the text); 100,000
                                     #   replications per cell (writes data/table_1.dat)
scripts/table_1_core.jl              # shared functions (v_I, omega, truncated exponential
                                     #   sampler); run alone, the asymptotic test with the
                                     #   first-order CLT predictions (data/table_1_clt.dat)
scripts/run_all.jl                   # all of the above
R/power_comparison.R                 # Section 6: comparison with the universal inference
                                     #   test of Dunn et al. (2,000 replications; about
                                     #   25 minutes; writes data/power_comparison.dat)
R/swiss_example.R                    # Section 6: data example on the R dataset `swiss`
                                     #   (writes data/swiss_example.dat)
figures/figure_N.pdf                 # the figures
data/                                # simulation output for Figures 4 and 5, Table 1 and
                                     #   Section 6
```

## Running

With Julia 1.10 or later, from the repository root:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia -t auto --project=. scripts/run_all.jl
```

The simulation for Figures 4 and 5 takes about 3 minutes with 8 threads,
Table 1 several minutes, and the other figures seconds.

The two R scripts need R 4 or later and the package
[`logcondens`](https://CRAN.R-project.org/package=logcondens):

```bash
Rscript -e 'install.packages("logcondens")'
Rscript R/swiss_example.R
Rscript R/power_comparison.R
``` Each script can also be run on its own,
for example `julia --project=. scripts/figure_2.jl`. The plotting script
for Figures 4 and 5 reads the data written by the simulation script.
