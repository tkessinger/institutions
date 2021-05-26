#!/usr/bin/env julia

## plt_fixation_basins.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Plot fixation probability from Institutions simulations.
## This script corresponds to supplementary figures 1 and 2.

using CSV, PyPlot, Statistics

# load simulation output as a dataframe
runs = CSV.read("output/fixation_paper_basins.csv")

# get unique values from the runs dataframe
N = sort(unique(runs[:N]))[1]
num_samples = sort(unique(runs[:num_trials]))[1]

b_vals = unique(runs[:b])
Q_vals = unique(runs[:Q])
q_vals = unique(runs[:q])
E_vals = unique(runs[:E])
norms = ["stern judging", "simple standing", "scoring", "shunning"]

for (bi, b) in enumerate(b_vals)
	for run_type in unique(runs[:run_type])
		for indep in [1, 0]

			# dicts to store fixation/extinction events
			successes = Dict{Tuple{String, Float64, Int64, Float64},Int64}()
			failures = Dict{Tuple{String, Float64, Int64, Float64},Int64}()

			param_combs = collect(Base.product(norms, E_vals, Q_vals, q_vals))

			for (pi, param_comb) in enumerate(param_combs)
			    norm, E, Q, q = param_comb
			    successes[param_comb] = 0
				failures[param_comb] = 0
				# filter for correct parameter values
			    tmp_runs = runs[(runs[:Q] .== Q) .& (runs[:b] .== b) .& (runs[:E] .== E) .& (runs[:q] .== q) .& (runs[:reputation_norm] .== norm) .& (runs[:independent_board] .== indep) .& (runs[:run_type] .== run_type),:]
			    for (ri, run) in enumerate(eachrow(tmp_runs[:,:]))
					successes[param_comb] += run[:successes]
					failures[param_comb] += run[:failures]
				end
			end

			param_combs = reshape(param_combs, length(param_combs))

			cs = ["tab:blue", "tab:green", "tab:red", "tab:orange"]

			sub_norms = ["stern judging", "simple standing", "scoring", "shunning"]

			fig, axs = plt.subplots(length(E_vals), length(sub_norms), figsize=(3*length(sub_norms), 3*length(E_vals)), sharey="row", sharex="col")
			for (Ei, E) in enumerate(E_vals)
				for (ni, norm) in enumerate(sub_norms)

					ax = axs[Ei, ni]
					title_string = "$norm"

					for (Qi, Q) in enumerate(Q_vals)
						fixations = [successes[norm, E, Q, q] for q in q_vals]
						extinctions = [failures[norm, E, Q, q] for q in q_vals]
						if Q == 2
							# hack to produce the "jump" in our Q = 2 plots
							mod_q_vals = vcat(q_vals[1:2], [0.499999], q_vals[3:end])
							fixations = vcat(fixations[1:2], [fixations[2]], fixations[3:end])
							extinctions = vcat(extinctions[1:2], [extinctions[2]], extinctions[3:end])
							p_val = 1.0*fixations./(fixations .+ extinctions)
							e_bars = sqrt.(p_val.*(1.0 .- p_val)./(fixations .+ extinctions))
							ax.errorbar(mod_q_vals, p_val, e_bars,
								label = "Q = $Q", c = cs[Qi])
						else
							p_val = 1.0*fixations./(fixations .+ extinctions)
							e_bars = sqrt.(p_val.*(1.0 .- p_val)./(fixations .+ extinctions))
							ax.errorbar(q_vals, p_val, 2*e_bars,
								label = "Q = $Q", c = cs[Qi])
						end
					end
					ax.hlines(0.5, 0, 1, linestyle = "--")

					if Ei == 1
						ax.set_title(title_string)
					end
					if ni == 1
						ax.set_ylabel("E = $E")
						if Ei == 1
							ax.legend(loc=2)
						end
					end
					ax.set_ylim([0,1])
					ax.set_xlim([0,1])

				end
			end
			titlestring = ""
			if indep == 1
				titlestring *= "external board"
			else
				titlestring *= "internal board"
			end
			titlestring *= ", b = $b"
			fig.suptitle("$titlestring")
			fig.text(0.5,0.04, "institution strictness " * L" q", ha="center", va="center", fontsize="16")
			fig.text(0.02,0.5, "mean fixation probability", ha="center", va="center", fontsize="16", rotation="90")
			fig.tight_layout(rect=[0.02, 0.04, 1, 0.96])
			display(fig)
			plt.savefig("figures/fig5_basins_board_$(indep)_var_b_$b.pdf")
		end
	end
end
