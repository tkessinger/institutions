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

strategies = [1,2,3]

verbose = false

w = 1.0

q = 0.5
Q = 2

E = 0.9
#rep_norm = "simple standing"
rep_norms = ["stern judging", "simple standing", "scoring", "shunning"]
rep_norm = rep_norms[1]

game = Game(b, c, w, u_s, u_p, u_a, "pc")

successful = false
trials = 0

while successful == false
        global trials += 1
        println("$trials")

        global strat_freqs = []
        global follower_freqs = []

        pop = fixation_population(N, E, Q, q, game, rep_norm, strategies, verbose)
        for i in 1:500
                evolve!(pop)
                push!(strat_freqs, get_freqs(pop))
                push!(follower_freqs, 1.0*sum(pop.is_follower)/pop.N)
        end
        if sum([pop.strategies[x] .== 3 for x in 1:pop.N]) > 0
                pop.is_follower[rand(filter(x -> pop.strategies[x] .== 3, 1:pop.N))] = 1
        else
                x = rand(1:pop.N)
                pop.strategies[x] = 3
                pop.is_follower[x] = 1
        end
        while sum(pop.is_follower) ∉ [0, 50]
        #for i in 1:1000
                evolve!(pop)
                push!(strat_freqs, get_freqs(pop))
                push!(follower_freqs, 1.0*sum(pop.is_follower)/pop.N)
        end
        if sum(pop.is_follower) == 50
                global successful = true
        end
end

strat_freqs = permutedims(hcat(strat_freqs...))
fig = plt.figure()
[plt.plot(strat_freqs[:,x], label="$x") for x in 1:4]
plt.plot(follower_freqs, ls = "--", label = "follower")
plt.ylim([0,1])
plt.tight_layout()
display(fig)
