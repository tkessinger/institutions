using Revise
using Institutions
using Statistics
using PyPlot

N = 4
Q = 1
q = 0.1

b = 5.0
c = 1.0

w = 1.0
u_s = 0.2
u_p = 0.0
u_a = 0.0

rep_norm = "scoring"

num_gens = 1

game = Game(b, c, w, u_s, u_p, u_a)

pop = Population(N, Q, q, game, rep_norm, [1,2,3,4], true)
for i in 1:num_gens
    evolve!(pop)
end
