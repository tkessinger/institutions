#!/usr/bin/env julia

## plt_test_results.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Plot results from Institutions simulations.
## Look at type frequency dynamics and cooperation
## as a function of c, M, and K.

using CSV, PyPlot, Statistics

# load simulation output as a dataframe
runs = CSV.read("output/test_institutions.csv")

# dicts to store frequencies
coop_freqs = Dict{Tuple{String, Int64, Float64},Array{Float64, 1}}()
strat_freqs = Dict{Tuple{String, Int64, Float64},Array{Array{Float64, 1}}}()


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
	#strat_freqs[param_comb]
    tmp_runs = runs[(runs[:Q] .== Q) .& (runs[:q] .== q) .& (runs[:reputation_norm] .== norm), :]
    for (ri, run) in enumerate(eachrow(tmp_runs[:,:]))
		println("$ri")
		push!(coop_freqs[param_comb], run[:coop_freqs])
        # coop_rows = split.(run[:total_cooperation][2:end-1], r", ")
		# rep_means_rows = split.(run[:reputation_means][2:end-1], r", ")
		# rep_stds_rows = split.(run[:reputation_stds][2:end-1], r", ")
		#set_rep_means_rows = split.(run[:set_reputation_means][2:end-1], ";")
    end
end

param_combs = reshape(param_combs, length(param_combs))
strat_ids = "compartmentalizer", "forgiving", "draconian"

plotcolors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple", "tab:brown"]

cs = ["blue", "green", "red"]

fig, axs = plt.subplots(1,4, figsize=(20,5), sharey="row")
for (ni, norm) in enumerate(norms)
	#fig = plt.figure()
	ax = axs[ni]
	title_string = "norm = $norm"

	coop_means = zeros(Float64, length(Q_vals), length(q_vals))
	coop_stds = zeros(Float64, length(Q_vals), length(q_vals))
	for (Qi, Q) in enumerate(Q_vals)
		for (qi, q) in enumerate(q_vals)
			coop_means[Qi, qi] = mean(coop_freqs[norm, Q, q])
			coop_stds[Qi, qi] = std(coop_freqs[norm, Q, q])
		end
		ax.errorbar(q_vals, coop_means[Qi,:],
			coop_stds[Qi,:]/sqrt(num_samples),
			label="Q = $Q",
			color=cs[Qi])
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
