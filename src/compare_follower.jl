using Revise
using Institutions
using Statistics
using PyPlot

N = 50

num_trials = 1

b = 5.0
c = 1.0

u_s = 0.02
#u_s = 0.025
u_p = 0.02
u_a = 0.02

strategies = [3]

verbose = false

w = 0.0

q = 0.1
Q = 2

E = 0.0
#rep_norm = "simple standing"
rep_norms = ["stern judging", "simple standing", "scoring", "shunning"]
rep_norm = rep_norms[4]

game = Game(b, c, w, u_s, u_p, u_a, "pc")

successful = false
trials = 0

follower_reputations = []
indep_reputations = []
follower_freqs = []

pop = fixation_population(N, E, Q, q, game, rep_norm, strategies, verbose)
pop.is_follower[rand(1:N,N÷2)] .= 1
for i in 1:2500
        if rand() < pop.game.u_s
                x = rand(1:N)
                pop.is_follower[x] = 1 - pop.is_follower[x]
        end
        evolve!(pop)
        push!(follower_freqs, 1.0*sum(pop.is_follower)/pop.N)
        push!(follower_reputations, sum(pop.priv_reputations[:, pop.is_follower .== 1])/(pop.N*sum(pop.is_follower)))
        push!(indep_reputations, sum(pop.priv_reputations[:, pop.is_follower .== 0])/(pop.N*sum(1.0 .- pop.is_follower)))
end

fig = plt.figure()
plt.plot(follower_reputations, label = "follower reputation")
plt.plot(indep_reputations, label = "independent reputation")
plt.plot(follower_freqs, ls = "--", label = "follower frequency")
plt.ylim([0,1])
plt.legend(loc = 4)
plt.tight_layout()
display(fig)
