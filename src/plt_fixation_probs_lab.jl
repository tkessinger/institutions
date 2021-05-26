#!/usr/bin/env julia

## plt_fixation_probs.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Plot fixation probability from Institutions simulations.

using CSV, PyPlot, Statistics

# load simulation output as a dataframe
runs = CSV.read("output/fixation_paper_independent_redo.csv")

# dicts to store frequencies
successes = Dict{Tuple{String, Float64, Int64, Float64},Int64}()
failures = Dict{Tuple{String, Float64, Int64, Float64},Int64}()
reputations = Dict{Tuple{String, Float64, Int64, Float64},Float64}()

# get unique values from the runs dataframe
N = sort(unique(runs[:N]))[1]
num_samples = sort(unique(runs[:num_trials]))[1]

Q_vals = sort(unique(runs[:Q]))
q_vals = sort(unique(runs[:q]))
E_vals = unique(runs[:E])
norms = ["stern judging", "simple standing", "scoring", "shunning"]

param_combs = collect(Base.product(norms, E_vals, Q_vals, q_vals))

for (pi, param_comb) in enumerate(param_combs)
    norm, E, Q, q = param_comb
    #println("$param_comb")
    successes[param_comb] = 0
	failures[param_comb] = 0
	reputations[param_comb] = 0.0
	#strat_freqs[param_comb]
    tmp_runs = runs[(runs[:Q] .== Q) .& (runs[:E] .== E) .& (runs[:q] .== q) .& (runs[:reputation_norm] .== norm), :]
	println("$param_comb, $(size(tmp_runs, 1))")
    for (ri, run) in enumerate(eachrow(tmp_runs[:,:]))
		#println("$ri")
		# push!(successes[param_comb], run[:successes])
		# push!(failures[param_comb], run[:failures])
		# push!(reputations[param_comb], run[:reputations])
		successes[param_comb] += run[:successes]
		failures[param_comb] += run[:failures]
		reputations[param_comb] += run[:reputations]/size(tmp_runs, 1)
	end
end

param_combs = reshape(param_combs, length(param_combs))

for (pi, param_comb) in enumerate(param_combs)
	println("$param_comb, $(successes[param_comb]), $(failures[param_comb]), $(reputations[param_comb])")
end

plotcolors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple", "tab:brown"]

cs = ["blue", "green", "red"]

sub_norms = ["stern judging", "simple standing", "scoring", "shunning"]

fig, axs = plt.subplots(length(E_vals), length(sub_norms), figsize=(5*length(E_vals), length(sub_norms)), sharey="row", sharex="col")
for (Ei, E) in enumerate(E_vals)
	for (ni, norm) in enumerate(sub_norms)
		#fig = plt.figure()
		ax = axs[Ei, ni]
		title_string = "E = $E, $norm"

		for (Qi, Q) in enumerate(Q_vals)
			fixations = [successes[norm, E, Q, q] for q in q_vals]
			extinctions = [failures[norm, E, Q, q] for q in q_vals]
			if Q == 2
				mod_q_vals = vcat(q_vals[1:2], [0.499999], q_vals[3:end])
				fixations = vcat(fixations[1:2], [fixations[2]], fixations[3:end])
				extinctions = vcat(extinctions[1:2], [extinctions[2]], extinctions[3:end])
				ax.plot(mod_q_vals, 1.0*fixations./(fixations .+ extinctions),
					label = "Q = $Q", c = cs[Qi])
			else
				ax.plot(q_vals, 1.0*fixations./(fixations .+ extinctions),
					label = "Q = $Q", c = cs[Qi])
			end
		end
		ax.hlines(1.0/N, 0, 1, linestyle = "--")
		equ_reputation = sum([sum([reputations[norm, E, Q, q] for q in q_vals]) for Q in Q_vals])/(length(q_vals)*length(Q_vals))
		ax.vlines(equ_reputation, 0, 1, linestyle = "--")
		ax.set_xlabel("q")
		if ni == 1 && Ei == 1
			ax.legend(loc=3)
			ax.set_ylabel("cooperation")
		end
		ax.set_ylim([0,1])
		ax.set_xlim([0,1])

		ax.set_title(title_string)
	end
end
fig.tight_layout(rect=[0, 0.03, 1, 0.96])
display(fig)
