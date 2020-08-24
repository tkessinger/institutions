#!/usr/bin/env julia

## plt_test_results.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Plot results from Institutions simulations.
## Look at cooperation frequency
## as a function of Q, q, and social norm.

using CSV, PyPlot, Statistics

# load simulation output as a dataframe
runs = CSV.read("output/test_institutions_mu_fixed_cumulative.csv")

# dicts to store frequencies
coop_freqs = Dict{Tuple{String, Int64, Float64},Array{Float64, 1}}()
reps = Dict{Tuple{String, Int64, Float64},Array{Float64, 1}}()
strat_freqs = Dict{Tuple{String, Int64, Float64},Array{Array{Float64, 1}}}()

E_coop_freqs = Dict{Tuple{String, Float64},Array{Float64, 1}}()


empathy_runs = CSV.read("output/test_empathy_cumulative_RDisc.csv")

# get unique values from the runs dataframe
N = sort(unique(runs[:N]))[1]

Q_vals = unique(runs[:Q])
q_vals = unique(runs[:q])
#norms = unique(runs[:reputation_norm])
norms = ["stern judging", "simple standing", "scoring", "shunning"]

num_samples = sort(unique(runs[:num_trials]))[1]

param_combs = collect(Base.product(norms, Q_vals, q_vals))

for (pi, param_comb) in enumerate(param_combs)
    norm, Q, q = param_comb
	println("$param_comb")
    #println("$param_comb")
    coop_freqs[param_comb] = Float64[]
	reps[param_comb] = Float64[]
	strat_freqs[param_comb] = Array{Float64,1}[]
	#strat_freqs[param_comb]
    tmp_runs = runs[(runs[:Q] .== Q) .& (runs[:q] .== q) .& (runs[:reputation_norm] .== norm), :]
    for (ri, run) in enumerate(eachrow(tmp_runs[:,:]))
		println("$ri")
		push!(coop_freqs[param_comb], run[:coop_freqs])
		push!(reps[param_comb], run[:reputations])
        freqs_rows = split.(run[:strat_freqs][2:end-1], r", ")
		push!(strat_freqs[param_comb], map(x -> parse(Float64,x), freqs_rows))
		# rep_means_rows = split.(run[:reputation_means][2:end-1], r", ")
		# rep_stds_rows = split.(run[:reputation_stds][2:end-1], r", ")
		#set_rep_means_rows = split.(run[:set_reputation_means][2:end-1], ";")
    end
end

E_vals = sort(unique(empathy_runs[:E]))
empathy_samples = unique(empathy_runs[:num_trials])

E_params = collect(Base.product(norms, E_vals))

for (pi, param_comb) in enumerate(E_params)
    norm, E = param_comb
	println("$param_comb")
    #println("$param_comb")
    E_coop_freqs[param_comb] = Float64[]
    tmp_runs = empathy_runs[(empathy_runs[:E] .== E) .& (empathy_runs[:reputation_norm] .== norm), :]
    for (ri, run) in enumerate(eachrow(tmp_runs[:,:]))
		println("$ri")
		push!(E_coop_freqs[param_comb], run[:coop_freqs])
    end
end

param_combs = reshape(param_combs, length(param_combs))
E_param_combs = reshape(E_params, length(E_params))

plotcolors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple", "tab:brown"]

cs = ["blue", "green", "red"]

e_cs = ["yellow", "black"]

fig, axs = plt.subplots(1,4, figsize=(20,5), sharey="row")
for (ni, norm) in enumerate(norms)
	#fig = plt.figure()
	ax = axs[ni]
	title_string = "norm = $norm"

	E_coop_means = zeros(Float64, length(E_vals))

	coop_means = zeros(Float64, length(Q_vals), length(q_vals))
	coop_stds = zeros(Float64, length(Q_vals), length(q_vals))

	rep_means = zeros(Float64, length(Q_vals), length(q_vals))
	rep_stds = zeros(Float64, length(Q_vals), length(q_vals))

	for (Qi, Q) in enumerate(Q_vals)
		for (qi, q) in enumerate(q_vals)
			coop_means[Qi, qi] = mean(coop_freqs[norm, Q, q])
			coop_stds[Qi, qi] = std(coop_freqs[norm, Q, q])
			rep_means[Qi, qi] = mean(reps[norm, Q, q])
			rep_stds[Qi, qi] = std(reps[norm, Q, q])
		end
		ax.errorbar(q_vals, coop_means[Qi,:],
			coop_stds[Qi,:]/sqrt(num_samples),
			label="Q = $Q",
			color=cs[Qi])
		#ax.errorbar(q_vals, rep_means[Qi,:],
		#	rep_stds[Qi,:]/sqrt(num_samples),
			#label="Q = $Q",
		#	ls = "--",
		#	color=cs[Qi])
	end
	for (Ei, E) in enumerate(E_vals)
		E_coop_means[Ei] = mean(E_coop_freqs[norm, E])
		ax.hlines(E_coop_means[Ei], 0, 1, linestyle="--", label = "empathy = $E", color = e_cs[Ei])
	end
	ax.set_xlabel("q")
	if ni == 1
		ax.legend(loc=3)
		ax.set_ylabel("cooperation")
	end
	ax.set_ylim([0,1])
	ax.set_xlim([0,1])

	ax.set_title(title_string)
end
fig.tight_layout(rect=[0, 0.03, 1, 0.96])
display(fig)

# fig, axs = plt.subplots(length(Q_vals), length(norms), figsize=(20,20), sharey="row", sharex="col")
# for (Qi, Q) in enumerate(Q_vals)
# 	for (ni, norm) in enumerate(norms)
# 		#fig = plt.figure()
# 		ax = axs[Qi, ni]
# 		#title_string = "norm = $norm"
#
# 		strat_means = zeros(Float64, 4, length(q_vals))
# 		strat_stds = zeros(Float64, 4, length(q_vals))
#
# 		for (qi, q) in enumerate(q_vals)
# 			strat_means[:, qi] = [mean(hcat(strat_freqs[norm, Q, q]...)[x,:]) for x in 1:4]
# 			strat_stds[:, qi] = [std(hcat(strat_freqs[norm, Q, q]...)[x,:]) for x in 1:4]
# 		end
# 		for i in 1:4
# 			ax.errorbar(q_vals, strat_means[i,:],
# 				strat_stds[i,:]/sqrt(num_samples),
# 				label="$i")
# 		end
# 		#color=cs[Qi])
# 		if Qi == 1
# 			ax.set_title("$norm")
# 		end
# 		if Qi == length(Q_vals)
# 			ax.set_xlabel("q")
# 		end
# 		if ni == 1
# 			ax.legend(loc=3)
# 			ax.set_ylabel("Q = $Q")
# 		end
# 		ax.set_ylim([0,1])
# 		ax.set_xlim([0,1])
#
# 		#ax.set_title(title_string)
# 	end
# end
# fig.tight_layout(rect=[0, 0.03, 1, 0.96])
# display(fig)
