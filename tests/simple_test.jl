#!/usr/bin/env julia

## simple_test.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Test script for Institutions.
## Simulate a small population and plot
## strategy frequencies and reputations
## for stern judging and shunning.

using Revise
using Institutions
using Statistics
using PyPlot

N = 50 # population size
Q = 2 # institution size
q = 0.7 # strictness parameter

b = 5.0 # benefit to cooperation
c = 1.0 # cost to cooperation

w = 1.0 # strength of selection
u_s = 0.02 # mutation rate
u_p = 0.02 # performance error (e1)
u_a = 0.02 # assessment error (e2)

permitted_strategies = [1,2,3] # AllC, AllD, Disc

num_gens = 10000

for rep_norm in ["stern judging", "shunning"]

    # initialize game object
    game = Game(b, c, w, u_s, u_p, u_a)

    # array to store strategy frequencies
    strategies = Array{Float64,1}[]
    reputation = Float64[]

    pop = institution_population(N, Q, q, game, rep_norm, permitted_strategies, true)
    for i in 1:num_gens
        # evolve population
        evolve!(pop)
        push!(strategies, get_freqs(pop))
        push!(reputation, get_pub_reputations(pop))
    end

    # reshape strategies
    strategies = permutedims(hcat(strategies...))
    strat_label = ["AllC", "AllD", "Disc"]

    # plot frequencies and reputation
    fig = plt.figure()
    for i in 1:3
        plt.plot(collect(1:num_gens), strategies[:,i], label = strat_label[i])
    end
    plt.plot(collect(1:num_gens), reputation, ls = "--", label = "reputation", c = "k")
    plt.legend(loc=2)
    plt.ylim([0,1])
    plt.xlabel("time")
    plt.ylabel("frequency")
    plt.title("$rep_norm")
    plt.tight_layout()
    display(fig)
end

# with this parameter combination, you should see:
# 1. reputations are high under stern judging
# UNLESS ALLD is dominating. then reputations fluctuate.
# 2. reputations are low under shunning no matter what.
