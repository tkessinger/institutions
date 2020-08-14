using Revise
using Institutions
using Statistics
using PyPlot

N = 50
Q = 50

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
    #for rep_norm in ["scoring"]
        fig, axs = plt.subplots(5,1, figsize=(10,20), sharex="col")
        for (qi, q) in enumerate([0.1, 0.3, 0.5, 0.7, 0.9])
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

            ax = axs[qi]
            ax.plot(coop, ls = "--", c = "black", label="coop")
            ax.plot(reputations, ls = "-", c = "yellow", label="rep")
            [ax.plot(strat_freqs[:,x], label="$x") for x in 1:4]
            ax.set_ylim([0,1])
            ax.set_xlim([0,num_gens])
            #ax.title("norm = $rep_norm, w = $w")
            ax.set_ylabel("q = $q")
            if qi == 1
                ax.legend(loc=2)
                ax.set_title("norm = $rep_norm, Q = $Q, w = $w")
            end
            if qi == 5
                ax.set_xlabel("time")
            end
        end
        plt.tight_layout()
        display(fig)
        plt.savefig("multi_norm_$(replace(rep_norm, " " => "_"))_w_$(w)_Q_$(Q).pdf")
    end
end