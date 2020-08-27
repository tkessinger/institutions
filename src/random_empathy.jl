using Revise
using Empathy
using Statistics
using PyPlot
using Random

Random.seed!(1)

N = 3
E = 1.0

b = 5.0
c = 1.0

w = 1.0
u_s = 0.2
u_p = 0.1
u_a = 0.1

rep_norm = "shunning"

num_gens = 1

game = Game(b, c, w, u_s, u_p, u_a)

pop = Population(N, E, game, rep_norm, [1,2,3,4], true)

println("initial strategies are $(pop.strategies)")
println("initial private reputations are $(pop.priv_reputations)")
for i in 1:num_gens
    evolve!(pop)
end
println("final strategies are $(pop.strategies)")
println("final private reputations are $(pop.priv_reputations)")
