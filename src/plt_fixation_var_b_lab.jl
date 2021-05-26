#!/usr/bin/env julia

## plt_fixation_probs.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Plot fixation probability from Institutions simulations.

using CSV, PyPlot, Statistics

# load simulation output as a dataframe
runs = CSV.read("output/fixation_paper_var_b.csv")

# get unique values from the runs dataframe
N = sort(unique(runs[:N]))[1]
num_samples = sort(unique(runs[:num_trials]))[1]

b_vals = unique(runs[:b])
Q_vals = unique(runs[:Q])
q_vals = unique(runs[:q])
E_vals = unique(runs[:E])
norms = ["stern judging", "simple standing", "scoring", "shunning"]

fig, axs = plt.subplots(2, length(norms), figsize=(10, length(sub_norms)), sharey="row", sharex="col")

runs = CSV.read("output/fixation_paper_var_b.csv")
# dicts to store frequencies
successes = Dict{Tuple{String, Float64, Int64, Float64},Int64}()
failures = Dict{Tuple{String, Float64, Int64, Float64},Int64}()
#reputations = Dict{Tuple{String, Float64, Int64, Float64},Float64}()


param_combs = collect(Base.product(norms, [E_vals[1]], Q_vals, q_vals))

b = 60.0
indep = true
run_type = "institution"

for (pi, param_comb) in enumerate(param_combs)
    norm, E, Q, q = param_comb
    #println("$param_comb")
    successes[param_comb] = 0
	failures[param_comb] = 0
	#reputations[param_comb] = 0.0
	#strat_freqs[param_comb]
    tmp_runs = runs[(runs[:Q] .== Q) .& (runs[:b] .== b) .& (runs[:E] .== E) .& (runs[:q] .== q) .& (runs[:reputation_norm] .== norm) .& (runs[:independent_board] .== indep) .& (runs[:run_type] .== run_type),:]
	#println("$param_comb, $(size(tmp_runs, 1))")
    for (ri, run) in enumerate(eachrow(tmp_runs[:,:]))
		#println("$ri")
		# push!(successes[param_comb], run[:successes])
		# push!(failures[param_comb], run[:failures])
		# push!(reputations[param_comb], run[:reputations])
		#println("$param_comb")
		successes[param_comb] += run[:successes]
		failures[param_comb] += run[:failures]
		#reputations[param_comb] += run[:reputations]/size(tmp_runs, 1)
	end
end

param_combs = reshape(param_combs, length(param_combs))

# for (pi, param_comb) in enumerate(param_combs)
# 	println("$param_comb")
# 	println("$(successes[param_comb]), $(failures[param_comb])")
# end

plotcolors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple", "tab:brown"]

cs = ["blue", "green", "red"]

sub_norms = ["stern judging", "simple standing", "scoring", "shunning"]

b = 60.0
bi = 2
indep = 1

for (Ei, E) in enumerate([E_vals[1]])
	for (ni, norm) in enumerate(sub_norms)
		#fig = plt.figure()
		ax = axs[2, ni]
		title_string = "b = $b, $norm"

		for (Qi, Q) in enumerate(Q_vals)
			fixations = [successes[norm, E, Q, q] for q in q_vals]
			extinctions = [failures[norm, E, Q, q] for q in q_vals]
			if Q == 2
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
		#equ_reputation = sum([sum([reputations[norm, E, Q, q] for q in q_vals]) for Q in Q_vals])/(length(q_vals)*length(Q_vals))
		#ax.vlines(equ_reputation, 0, 1, linestyle = "--")
		ax.set_xlabel("q")
		if ni == 1
			ax.set_ylabel("fixation probability")
			if bi == 1
				ax.legend(loc=2)
			end
		end
		ax.set_ylim([0,1])
		ax.set_xlim([0,1])

		ax.set_title(title_string)
	end
end

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

b = 5.0
bi = 1
indep = 1

for (Ei, E) in enumerate([E_vals[1]])
	for (ni, norm) in enumerate(sub_norms)
		#fig = plt.figure()
		ax = axs[1, ni]
		title_string = "b = $b, $norm"

		for (Qi, Q) in enumerate(Q_vals)
			fixations = [successes[norm, E, Q, q] for q in q_vals]
			extinctions = [failures[norm, E, Q, q] for q in q_vals]
			if Q == 2
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
		#equ_reputation = sum([sum([reputations[norm, E, Q, q] for q in q_vals]) for Q in Q_vals])/(length(q_vals)*length(Q_vals))
		#ax.vlines(equ_reputation, 0, 1, linestyle = "--")
		ax.set_xlabel("q")
		if ni == 1
			ax.set_ylabel("fixation probability")
			if bi == 1
				ax.legend(loc=2)
			end
		end
		ax.set_ylim([0,1])
		ax.set_xlim([0,1])

		ax.set_title(title_string)
	end
end
plt.tight_layout()

display(fig)
