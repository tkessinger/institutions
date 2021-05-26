using Revise
using Empathy
using Statistics
using PyPlot

N = 50
E = 0.0

b = 5.0
c = 1.0

w = 1.0
u_s = 0.025
u_p = 0.02
u_a = 0.02

strategies = [1,2,3]

rep_norm = "stern judging"

num_gens = 100000
verbose = false

for E in [0.0]
    fig, axs = plt.subplots(4,1, figsize = (10,10), sharex="col")
    for (ri,rep_norm) in enumerate(["stern judging", "simple standing", "scoring", "shunning"])

        game = Game(b, c, w, u_s, u_p, u_a, "pc")

        coop = Float64[]
        reputations = Float64[]
        strat_freqs = Array{Float64, 1}[]

        pop = Population(N, E, game, rep_norm, strategies, verbose)

        for i in 1:num_gens
            evolve!(pop)
            push!(coop, sum(pop.prev_actions)/(N*(N-1)))
            push!(reputations, get_reputations(pop))
            push!(strat_freqs, get_freqs(pop))
        end

        strat_freqs = permutedims(hcat(strat_freqs...))

        ax = axs[ri]

        ax.plot(coop, ls = "--", c = "black", label="coop")
        ax.plot(reputations, ls = "-", c = "yellow", label="rep")
        [ax.plot(strat_freqs[:,x], label="$x") for x in 1:4]
        ax.set_ylim([0,1])
        ax.set_xlim([0,num_gens])
        ax.set_title("norm = $rep_norm, E = $E, w = $w")
        if ri == 1
            ax.legend(loc=2)
        end
        if ri == 4
            ax.set_xlabel("time")
        end
        ax.set_ylabel("frequency")
    end
    plt.tight_layout()
    display(fig)
    plt.savefig("empathy_sims_E_$(E)_2.pdf")
end
