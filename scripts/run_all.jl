#=
run_all.jl

Reproduces all figures. Run from the repository root as

    julia -t auto --project=. scripts/run_all.jl

The simulation for Figures 4 and 5 takes about 3 minutes with 8 threads;
the other figures take seconds.
=#

for script in ("figures_1_and_6.jl", "figure_2.jl", "figure_3.jl",
               "figures_4_and_5_simulate.jl", "figures_4_and_5_plot.jl")
    println("== ", script)
    include(joinpath(@__DIR__, script))
end
