using Revise
using Institutions
using Statistics
using PyPlot

N = 50

num_trials = 1

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
#rep_norm = "simple standing"
rep_norms = ["stern judging", "simple standing", "scoring", "shunning"]

game = Game(b, c, w, u_s, u_p, u_a, "pc")

for rep_norm in rep_norms
        for i in 1:num_trials
                equ_reputations = []
                pop = fixation_population(N, E, Q, q, game, rep_norm, strategies, verbose)
                for i in 1:1000
                        evolve!(pop)
                        push!(equ_reputations, sum(pop.priv_reputations)/N^2)
                end
                fig = plt.figure()
                plt.plot(equ_reputations)
                plt.ylim([0,1])
                plt.title("$rep_norm, E = $E")
                plt.tight_layout()
                display(fig)
        end
end
