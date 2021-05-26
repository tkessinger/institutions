using Revise
using Institutions
using Statistics
using PyPlot

N = 25
#Q = 50

b = 5.0
c = 1.0

u_s = 0.025
u_p = 0.04
u_a = 0.04

strategies = [1,2,3,4]

verbose = false

rep_norm = "shunning"
w = 1.0

num_gens = 200000
Q = 1
q = 0.1

game = Game(b, c, w, u_s, u_p, u_a, "pc")

coop = Float64[]
reputations = Float64[]
strat_reps = Array{Float64, 1}[]
strat_freqs = Array{Float64, 1}[]

pop = empathy_population(N, 1.0, game, rep_norm, strategies, true)
for i in 1:num_gens
    evolve!(pop)
    push!(coop, sum(pop.prev_actions)/N^2)
    push!(reputations, get_priv_reputations(pop))
    push!(strat_freqs, get_freqs(pop))
    push!(strat_reps, get_strat_priv_reputations(pop))
end

strat_freqs = permutedims(hcat(strat_freqs...))
strat_reps = permutedims(hcat(strat_reps...))


fig, axs = plt.subplots(2,2,figsize=(10,10),sharex="col",sharey="row")

cs = ["blue", "orange", "green", "red"]

println("E: coop = $(mean(coop)), $(std(coop)), rep = $(mean(reputations)), $(std(reputations))")

ax = axs[1]

ax.plot(coop, ls = "--", c = "black", label="coop")
ax.plot(reputations, ls = "-", c = "yellow", label="rep")
[ax.plot(strat_freqs[:,x], label="$x", c = cs[x]) for x in 1:4]
#[plt.plot(strat_reps[:,x], c = cs[x], ls = "--") for x in 1:4]
ax.set_ylim([0,1])
ax.set_xlim([0,num_gens])
#ax.title("norm = $rep_norm, w = $w")
ax.legend(loc=2)
ax.set_title("empathy")
ax.set_ylabel("frequency")


ax = axs[2]

# ax.plot(coop, ls = "--", c = "black", label="coop")
# ax.plot(reputations, ls = "-", c = "yellow", label="rep")
[ax.plot(strat_reps[:,x], label="$x", c = cs[x], ls = "--") for x in 1:4]
ax.set_ylim([0,1])
ax.set_xlim([0,num_gens])
#ax.title("norm = $rep_norm, w = $w")
ax.legend(loc=2)
ax.set_ylabel("reputation")
ax.set_xlabel("time")

coop = Float64[]
reputations = Float64[]
strat_freqs = Array{Float64, 1}[]
strat_reps = Array{Float64, 1}[]

pop = institution_population(N, Q, q, game, rep_norm, strategies, true)
for i in 1:num_gens
    evolve!(pop)
    push!(coop, sum(pop.prev_actions)/N^2)
    push!(reputations, get_pub_reputations(pop))
    push!(strat_freqs, get_freqs(pop))
    push!(strat_reps, get_strat_pub_reputations(pop))
end

strat_freqs = permutedims(hcat(strat_freqs...))
strat_reps = permutedims(hcat(strat_reps...))

println("I: coop = $(mean(coop)), $(std(coop)), rep = $(mean(reputations)), $(std(reputations))")

ax =axs[3]

ax.plot(coop, ls = "--", c = "black", label="coop")
ax.plot(reputations, ls = "-", c = "yellow", label="rep")
[ax.plot(strat_freqs[:,x], label="$x", c = cs[x]) for x in 1:4]
#[plt.plot(strat_reps[:,x], c = cs[x], ls = "--") for x in 1:4]
ax.set_ylim([0,1])
ax.set_xlim([0,num_gens])
#ax.title("norm = $rep_norm, w = $w")
ax.set_title("institution")

ax = axs[4]

# ax.plot(coop, ls = "--", c = "black", label="coop")
# ax.plot(reputations, ls = "-", c = "yellow", label="rep")
[ax.plot(strat_reps[:,x], label="$x", c = cs[x], ls = "--") for x in 1:4]
ax.set_ylim([0,1])
ax.set_xlim([0,num_gens])
ax.set_xlabel("time")

#ax.title("norm = $rep_norm, w = $w")

plt.tight_layout()
plt.suptitle("$rep_norm",fontsize=18)
fig.subplots_adjust(top=0.95)
display(fig)
