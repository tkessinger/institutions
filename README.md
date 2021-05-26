# Institutions
Adherence to public institutions that foster cooperation
Arunas L. Radzvilavicius*, Taylor A. Kessinger*, and Joshua B. Plotkin
(*equal contribution)

# Overview
All code for this project is written for Julia 1.4.0.
The main simulation code is in module/Institutions.jl (the path to this file must be added to the user's LOAD_PATH variable, e.g., by adding push!(LOAD_PATH, "/path/to/module" to ~/.julia/config/startup.jl).
The module defines Game and Population types which must be instantiated with a constructor (several are provided).
Most important population information is stored as array attributes of the Population type.
Once the game is defined and initial population parameters, such as initial strategy composition, are specified, evolution is as simple as passing the population to the evolve! function.
A variety of population summary statistics, such as strategy frequencies, can easily be obtained.

In the src directory, we provide a test_sim.jl file that produces some simple simulation output.
The bulk of the simulation code is designed to be run on a computing cluster and relies on .json files in the submit directory.

To generate figure 3, run src/parallel_fig3_redo.jl and src/parallel_fig3_redo_E.jl, which store cooperation rates for institutional and empathetic populations (respectively) in output/paper_institutions_noRDisc.csv and output/paper_empathy_noRDisc.csv.
The file plt/plt_coop_freqs_discfreq.jl will plot the resulting data.

To generate figure 4 and supplementary figures 7-9, run src/parallel_fixation_variable_b.jl (output saved as output/fixation_paper_variable_b.csv) followed by plt/plt_fixation_var_b.jl.

To generate supplementary figure 1, run src/parallel_fixation_basins.jl (output saved as output/fixation_paper_basins.csv) followed by plt/plt_fixation_basins.jl.

To generate supplementary figure 14, run src/parallel_fixation_adherents.jl (output saved as output/fixation_paper_adherents.csv) followed by plt/plt_fixation_adherents.jl.

Finally, to generate supplementary figures 3-6, run src/parallel_fixation_multitype.jl (output saved as output/fixation_paper_multitype_final.jl) followed by plt/plt_fixation_multitype.jl.


