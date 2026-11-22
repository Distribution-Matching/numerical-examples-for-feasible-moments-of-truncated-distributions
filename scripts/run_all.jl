#=
run_all.jl

Reproduces all figures, Table 1 and the numbers quoted in Sections 5 and 7.
Run from the repository root as

    julia -t auto --project=. scripts/run_all.jl

The simulation for Figure 4 takes about 3 minutes with 8 threads, Table 1
several minutes, and the other scripts seconds. The R scripts in R/
(comparison with the universal inference test, data example) and the
half-line size check (scripts/section_6_halfline_size.jl) are run
separately.
=#

for script in ("figures_1_and_5.jl", "figure_2.jl", "figure_3.jl",
               "figure_4_simulate.jl", "figure_4_plot.jl", "table_1.jl",
               "section_5_variance_peaks.jl", "section_5_inner_boundary.jl",
               "section_7_ifr_ceiling.jl")
    println("== ", script)
    # Each script runs in its own module, so that their global names do not
    # clash; `include` is defined there for scripts that include helpers.
    m = Module(Symbol(script))
    Core.eval(m, :(include(p) = Base.include($m, p)))
    Base.include(m, joinpath(@__DIR__, script))
end
