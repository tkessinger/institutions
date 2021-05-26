#!/usr/bin/env julia

## plt_fixation_multitype.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Plot fixation probability from Institutions simulations.
## This script corresponds to supplementary figures 3-6.

using CSV, PyPlot, Statistics

# load simulation output as a dataframe
runs = CSV.read("output/fixation_paper_multitype_final.csv")

# get unique values from the runs dataframe
N = sort(unique(runs[:N]))[1]
num_samples = sort(unique(runs[:num_trials]))[1]

Q_vals = unique(runs[:Q])
q_vals = unique(runs[:q])
E_vals = unique(runs[:E])
norms = ["stern judging", "simple standing", "scoring", "shunning"]

for run_type in unique(runs[:run_type])
	for indep in [1, 0]

		# dicts to store frequencies
		successes = Dict{Tuple{String, Float64, Int64, Float64},Int64}()
		failures = Dict{Tuple{String, Float64, Int64, Float64},Int64}()
		reputations = Dict{Tuple{String, Float64, Int64, Float64},Float64}()

		param_combs = collect(Base.product(norms, E_vals, Q_vals, q_vals))

		for (pi, param_comb) in enumerate(param_combs)
		    norm, E, Q, q = param_comb
		    successes[param_comb] = 0
			failures[param_comb] = 0
			reputations[param_comb] = 0.0
			# filter for correct parameter values
		    tmp_runs = runs[(runs[:Q] .== Q) .& (runs[:E] .== E) .& (runs[:q] .== q) .& (runs[:reputation_norm] .== norm) .& (runs[:independent_board] .== indep) .& (runs[:run_type] .== run_type),:]
		    for (ri, run) in enumerate(eachrow(tmp_runs[:,:]))
				successes[param_comb] += run[:successes]
				failures[param_comb] += run[:failures]
				reputations[param_comb] += run[:reputations]/size(tmp_runs, 1)
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
						# hack for the "jump" in our Q = 2 plots
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
				ax.hlines(1.0/N, 0, 1, linestyle = "--")
				equ_reputation = sum([sum([reputations[norm, E, Q, q] for q in q_vals]) for Q in Q_vals])/(length(q_vals)*length(Q_vals))
				ax.vlines(equ_reputation, 0, 1, linestyle = "--")
				#ax.set_xlabel("q")
				if ni == 1
					ax.set_ylabel("E = $E")
					if Ei == 1
						ax.legend(loc=2)
					end
				end
				if run_type == "AllD"
					ax.set_ylim([0,0.2])
					ax.set_yticks(collect(0:0.02:0.2))
				elseif run_type == "equal" || run_type == "AllC"
					ax.set_ylim([0,0.6])
				elseif run_type == "Disc"
					ax.set_ylim([0,0.4])
				end
				ax.set_xlim([0,1])
				if Ei == 1
					ax.set_title(title_string)
				end
			end
		end
		titlestring = ""
		if indep == 1
			titlestring *= "external board"
		else
			titlestring *= "internal board"
		end
		if run_type == "equal"
			titlestring *= ", central"
		else
			titlestring *= ", $run_type corner"
		end
		fig.suptitle("$titlestring")
		fig.text(0.5,0.04, "institution strictness q", ha="center", va="center", fontsize="16")
		fig.text(0.02,0.5, "fixation probability", ha="center", va="center", fontsize="16", rotation="90")
		fig.tight_layout(rect=[0.02, 0.04, 1, 0.96])
		display(fig)
		plt.savefig("figures/fig4_variant_board_$(indep)_initial_$run_type.pdf")
	end
end
