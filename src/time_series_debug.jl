using Revise
using Institutions
using Statistics
using PyPlot

N = 50
Q = 50
q = 0.3

b = 5.0
c = 1.0

u_s = 0.025
u_p = 0.02
u_a = 0.02

strategies = [1,2,3]

verbose = false

#rep_norm = "stern judging"
#w = 0.1

num_gens = 100000

for w in [0.1]
    for rep_norm in ["scoring", "shunning", "stern judging", "simple standing"]
        game = Game(b, c, w, u_s, u_p, u_a, "pc")

        coop = Float64[]
        reputations = Float64[]
        strat_freqs = Array{Float64, 1}[]

        pop = Population(N, Q, q, game, rep_norm, strategies, verbose)
        for i in 1:num_gens
            evolve!(pop)
            push!(coop, sum(pop.prev_actions)/(N*(N-1)))
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
    end
end
