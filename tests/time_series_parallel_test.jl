using Revise
using InstitutionsNew
using Statistics
using PyPlot
using Time

num_trials = 100

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

for j in 1:num_trials

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
end

for j in 1:num_trials

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
end
