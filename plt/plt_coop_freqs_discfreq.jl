#!/usr/bin/env julia

## plt_coop_freqs_discfreq.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Plot results from Institutions simulations.
## Look at cooperation frequency
## as a function of Q, q, and social norm.
## This script corresponds to figure 3.

using CSV, PyPlot, Statistics

# load simulation output as a dataframe
runs = CSV.read("output/paper_institutions_noRDisc.csv")

empathy_runs = CSV.read("output/paper_empathy_noRDisc.csv")

# dicts to store cooperation, frequencies, etc.
coop_freqs = Dict{Tuple{String, Int64, Float64},Array{Float64, 1}}()
reps = Dict{Tuple{String, Int64, Float64},Array{Float64, 1}}()
strat_freqs = Dict{Tuple{String, Int64, Float64},Array{Array{Float64, 1}}}()

E_coop_freqs = Dict{Tuple{String, Float64},Array{Float64, 1}}()
E_strat_freqs = Dict{Tuple{String, Float64},Array{Array{Float64, 1}}}()

# get unique values from the runs dataframe
N = sort(unique(runs[:N]))[1]

Q_vals = unique(runs[:Q])
q_vals = unique(runs[:q])
norms = ["stern judging", "simple standing", "scoring", "shunning"]

num_samples = sort(unique(runs[:num_trials]))[1]

param_combs = collect(Base.product(norms, Q_vals, q_vals))

for (pi, param_comb) in enumerate(param_combs)
    norm, Q, q = param_comb
    coop_freqs[param_comb] = Float64[]
	reps[param_comb] = Float64[]
	strat_freqs[param_comb] = Array{Float64,1}[]
	# filter for correct parameter values
    tmp_runs = runs[(runs[:Q] .== Q) .& (runs[:q] .== q) .& (runs[:reputation_norm] .== norm), :]
	println("$param_comb, $(size(tmp_runs, 1))")
    for (ri, run) in enumerate(eachrow(tmp_runs[:,:]))
		push!(coop_freqs[param_comb], run[:coop_freqs])
		push!(reps[param_comb], run[:reputations])
		# parse frequencies string into array
        freqs_rows = split.(run[:strat_freqs][2:end-1], r", ")
		push!(strat_freqs[param_comb], map(x -> parse(Float64,x), freqs_rows))
    end
end

E_vals = sort(unique(empathy_runs[:E]))
empathy_samples = unique(empathy_runs[:num_trials])

E_params = collect(Base.product(norms, E_vals))

for (pi, param_comb) in enumerate(E_params)
    norm, E = param_comb
    E_coop_freqs[param_comb] = Float64[]
	E_strat_freqs[param_comb] = Array{Float64,1}[]
    tmp_runs = empathy_runs[(empathy_runs[:E] .== E) .& (empathy_runs[:reputation_norm] .== norm), :]
	println("$param_comb, $(size(tmp_runs, 1))")
    for (ri, run) in enumerate(eachrow(tmp_runs[:,:]))
		push!(E_coop_freqs[param_comb], run[:coop_freqs])
		# parse frequencies string into array
		freqs_rows = split.(run[:strat_freqs][2:end-1], r", ")
		push!(E_strat_freqs[param_comb], map(x -> parse(Float64,x), freqs_rows))
    end
end

param_combs = reshape(param_combs, length(param_combs))
E_param_combs = reshape(E_params, length(E_params))

plotcolors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple", "tab:brown"]

cs = ["tab:blue", "tab:green", "tab:red", "tab:orange"]

e_cs = ["tab:purple", "black"]

sub_norms = ["stern judging", "simple standing", "scoring", "shunning"]

fig, axs = plt.subplots(1,length(sub_norms), figsize=(3*length(sub_norms),4), sharey="row")
for (ni, norm) in enumerate(sub_norms)
	ax = axs[ni]
	title_string = "$norm"

	E_coop_means = zeros(Float64, length(E_vals))

	coop_means = zeros(Float64, length(Q_vals), length(q_vals))
	coop_stds = zeros(Float64, length(Q_vals), length(q_vals))

	rep_means = zeros(Float64, length(Q_vals), length(q_vals))
	rep_stds = zeros(Float64, length(Q_vals), length(q_vals))

	strat_means = zeros(Float64, 4, length(q_vals))
	strat_stds = zeros(Float64, 4, length(q_vals))

	for (Qi, Q) in enumerate(Q_vals)
		for (qi, q) in enumerate(q_vals)
			coop_means[Qi, qi] = mean(coop_freqs[norm, Q, q])
			coop_stds[Qi, qi] = std(coop_freqs[norm, Q, q])
			rep_means[Qi, qi] = mean(reps[norm, Q, q])
			rep_stds[Qi, qi] = std(reps[norm, Q, q])

			strat_means[:, qi] = [mean(hcat(strat_freqs[norm, Q, q]...)[x,:]) for x in 1:4]
			strat_stds[:, qi] = [std(hcat(strat_freqs[norm, Q, q]...)[x,:]) for x in 1:4]
		end
		if Q == 2
			# hack to produce the "jump" visible in our Q = 2 plots
			mod_qvals = vcat(q_vals[1:2], [0.499999], q_vals[3:end])
			mod_means = vcat(coop_means[Qi,1:2], [coop_means[Qi,2]], coop_means[Qi,3:end])
			mod_stds = vcat(coop_stds[Qi,1:2], [coop_stds[Qi,2]], coop_stds[Qi,3:end])
			ax.errorbar(mod_qvals, mod_means,
				#2 for 99% confidence intervals
				2*mod_stds/sqrt(num_samples),
				label="Q = $Q",
				color=cs[Qi])
		else
			ax.errorbar(q_vals, coop_means[Qi,:],
				#2 for 99% confidence intervals
				2*coop_stds[Qi,:]/sqrt(num_samples),
				label="Q = $Q",
				color=cs[Qi])
		end
	end
	for (Ei, E) in enumerate(E_vals)
		E_coop_means[Ei] = mean(E_coop_freqs[norm, E])
		ax.hlines(E_coop_means[Ei], 0, 1, linestyle="--", label = "empathy = $E", color = e_cs[Ei])
	end
	if ni == 1
		ax.legend(loc=3)
		ax.set_ylabel("cooperation")
	end
	ax.set_ylim([0,1])
	ax.set_xlim([0,1])

	ax.set_title(title_string)
end
fig.text(0.5,0.04, "institution strictness q", ha="center", va="center", fontsize="16")
fig.tight_layout(rect=[0, 0.03, 1, 0.96])
display(fig)
plt.savefig("figures/fig3_paper_coop.pdf")

# uncomment the following to produce plots of
# strategy frequencies


# fig, axs = plt.subplots(length(Q_vals), length(sub_norms), figsize=(5*length(sub_norms),20), sharey="row", sharex="col")
# for (Qi, Q) in enumerate(Q_vals)
# 	for (ni, norm) in enumerate(sub_norms)
# 		ax = axs[Qi, ni]
# 		#title_string = "norm = $norm"
#
# 		coop_means = zeros(Float64, length(q_vals))
# 		coop_stds = zeros(Float64, length(q_vals))
#
# 		strat_means = zeros(Float64, 4, length(q_vals))
# 		strat_stds = zeros(Float64, 4, length(q_vals))
#
# 		for (qi, q) in enumerate(q_vals)
# 			strat_means[:, qi] = [mean(hcat(strat_freqs[norm, Q, q]...)[x,:]) for x in 1:4]
# 			strat_stds[:, qi] = [std(hcat(strat_freqs[norm, Q, q]...)[x,:]) for x in 1:4]
# 			coop_means[qi] = mean(coop_freqs[norm, Q, q])
# 			coop_stds[qi] = std(coop_freqs[norm, Q, q])
# 		end
# 		for i in 1:4
# 			ax.errorbar(q_vals, strat_means[i,:],
# 				strat_stds[i,:]/sqrt(num_samples),
# 				label="$i", c = cs[i])
# 			[ax.hlines(
# 				mean(hcat(E_strat_freqs[norm, E_val]...)[i,:]),
# 				minimum(q_vals), maximum(q_vals),
# 				linestyle = "--",
# 				color = cs[i])
# 				for E_val in [E_vals[1]]]
# 		end
# 		# ax.errorbar(q_vals, coop_means,
# 		# 	coop_stds/sqrt(num_samples),
# 		# 	label="coop",
# 		# 	color="teal",
# 		# 	ls = "--")
# 		#color=cs[Qi])
# 		if Qi == 1
# 			ax.set_title("$norm")
# 		end
# 		if Qi == length(Q_vals)
# 			ax.set_xlabel("institution strictness, q")
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
