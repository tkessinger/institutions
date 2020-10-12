using Revise
using Institutions
using Statistics
using PyPlot

N = 50

num_trials = 100

b = 5.0
c = 1.0

u_s = 0.0
#u_s = 0.025
u_p = 0.02
u_a = 0.02

strategies = [3]

verbose = false

w = 1.0

q = 0.5
Q = 1

E = 0.0
rep_norm = "stern judging"

game = Game(b, c, w, u_s, u_p, u_a, "pc")

global successes = 0
equilibrium_reputations = []

num_trials = 1000

for i in 1:num_trials
        #println("$i")
        pop = fixation_population(N, E, Q, q, game, rep_norm, strategies, verbose)
        for i in 1:100
                evolve!(pop)
        end
        init_reputation = sum(pop.priv_reputations)/N^2
        push!(equilibrium_reputations, init_reputation)
        pop.is_follower[rand(1:N)] = 1
        while sum(pop.is_follower) ∉ [0, pop.N]
                evolve!(pop)
                if sum(pop.is_follower) == pop.N
                        global successes += 1
                end
        end
        println("$i, $(pop.generation)")
end
