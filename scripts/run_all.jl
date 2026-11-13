#=
run_all.jl

Reproduces all figures and Table 1. Run from the repository root as

    julia -t auto --project=. scripts/run_all.jl

The simulation for Figures 4 and 5 takes about 3 minutes with 8 threads,
Table 1 about a minute, and the other figures seconds. The size check at
larger n quoted in the discussion of Table 1 (table_1_size_check.jl) takes
a few minutes and is not run here.
=#

for script in ("figures_1_and_6.jl", "figure_2.jl", "figure_3.jl",
               "figures_4_and_5_simulate.jl", "figures_4_and_5_plot.jl",
               "table_1.jl")
    println("== ", script)
    include(joinpath(@__DIR__, script))
end
