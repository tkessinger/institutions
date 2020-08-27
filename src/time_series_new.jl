using Revise
using InstitutionsNew
using Statistics
using PyPlot

N = 50
Q = 50
q = 0.5

b = 5.0
c = 1.0

u_s = 0.025
u_p = 0.02
u_a = 0.02

strategies = [1,2,3,4]

verbose = false

rep_norm = "shunning"
w = 1.0

num_gens = 10000

game = Game(b, c, w, u_s, u_p, u_a, "pc")

coop = Float64[]
reputations = Float64[]
strat_freqs = Array{Float64, 1}[]

pop = institution_population(N, Q, q, game, rep_norm, strategies, verbose)
#pop.strategies = 3*ones(Int64, pop.N)
for i in 1:num_gens
    evolve!(pop)
    push!(coop, sum(pop.prev_actions)/(N^2))
    push!(reputations, get_reputations(pop))
    push!(strat_freqs, get_freqs(pop))
end

strat_freqs = permutedims(hcat(strat_freqs...))

fig = plt.figure(figsize=(20,5))
plt.plot(coop, ls = "--", c = "black", label="coop")
plt.plot(reputations, ls = "-", c = "yellow", label="rep")
[plt.plot(strat_freqs[:,x], label="$x") for x in 1:4]
plt.ylim([0,1])
plt.xlim([0,num_gens])
plt.title("norm = $rep_norm, w = $w")
plt.legend(loc=2)
plt.tight_layout()
display(fig)
#plt.savefig("sample_norm_$(rep_norm)_w_$(w)_Q_$(Q)_q_$(q).pdf")


coop = Float64[]
reputations = Float64[]
strat_freqs = Array{Float64, 1}[]

E = 1.0

pop = empathy_population(N, E, game, rep_norm, strategies, verbose)
#pop.strategies = 3*ones(Int64, pop.N)
for i in 1:num_gens
    evolve!(pop)
    push!(coop, sum(pop.prev_actions)/(N^2))
    push!(reputations, get_reputations(pop))
    push!(strat_freqs, get_freqs(pop))
end

strat_freqs = permutedims(hcat(strat_freqs...))

fig = plt.figure(figsize=(20,5))
plt.plot(coop, ls = "--", c = "black", label="coop")
plt.plot(reputations, ls = "-", c = "yellow", label="rep")
[plt.plot(strat_freqs[:,x], label="$x") for x in 1:4]
plt.ylim([0,1])
plt.xlim([0,num_gens])
plt.title("norm = $rep_norm, w = $w")
plt.legend(loc=2)
plt.tight_layout()
display(fig)
#plt.savefig("sample_norm_$(rep_norm)_w_$(w)_Q_$(Q)_q_$(q).pdf")