using Revise
using Institutions
using Statistics
using PyPlot

N = 50
Q = 2
q = 0.9

b = 5.0
c = 1.0

w = 1.0
u_s = 0.025
u_p = 0.02
u_a = 0.02

verbose = false

rep_norm = "shunning"

num_gens = 10000

game = Game(b, c, w, u_s, u_p, u_a, "pc")

coop = Float64[]
reputations = Float64[]
strat_freqs = Array{Float64, 1}[]

pop = Population(N, Q, q, game, rep_norm, [1,2,3,4], verbose)
for i in 1:num_gens
    evolve!(pop)
    push!(coop, sum(pop.prev_actions)/(N*(N-1)))
    push!(reputations, get_reputations(pop))
    push!(strat_freqs, get_freqs(pop))
end

strat_freqs = permutedims(hcat(strat_freqs...))

fig = plt.figure()
plt.plot(coop, ls = "--", c = "black", label="coop")
plt.plot(reputations, ls = "-", c = "yellow", label="rep")
[plt.plot(strat_freqs[:,x], label="$x") for x in 1:4]
plt.ylim([0,1])
plt.legend(loc=2)
plt.tight_layout()
display(fig)
