using Revise
using InstitutionsNew
using Statistics
using PyPlot

N = 50

num_trials = 10

b = 5.0
c = 1.0

u_s = 0.025
u_p = 0.02
u_a = 0.02

strategies = [1,2,3,4]

verbose = false

w = 1.0

num_gens = 5000

q_vals = [0.1, 0.3, 0.5, 0.7, 0.9]

game = Game(b, c, w, u_s, u_p, u_a, "pc")

E_cs = ["yellow", "black"]

fig, axs = plt.subplots(1,4, figsize=(20, 5), sharey="row")

for (ri, rep_norm) in enumerate(["stern judging", "simple standing", "scoring", "shunning"])
    ax = axs[ri]
    E_coop_freqs = [[],[]]
    for (Ei, E) in enumerate([0.0, 1.0])
        for j in 1:num_trials
            coop = []
            pop = empathy_population(N, E, game, rep_norm, strategies, verbose)
            for i in 1:num_gens
                evolve!(pop)
                push!(coop, sum(pop.prev_actions/(N^2)))
            end
            push!(E_coop_freqs[Ei], mean(coop[num_gens÷2:end]))
        end
        ax.hlines(mean(E_coop_freqs[Ei]), 0, 1,
        linestyle="--", label="E = $E", color = E_cs[Ei])
    end

    for Q in [1, 2, 50]
        I_coop_freqs = [[],[],[],[],[]]
        for (qi, q) in enumerate(q_vals)
            println("$rep_norm, $Q, $q")
            for j in 1:num_trials
                coop = Float64[]

                pop = institution_population(N, Q, q, game, rep_norm, strategies, verbose)
                #pop.strategies = 3*ones(Int64, pop.N)
                for i in 1:num_gens
                    evolve!(pop)
                    push!(coop, sum(pop.prev_actions)/(N^2))
                end
                push!(I_coop_freqs[qi], mean(coop[num_gens÷2:end]))
            end
        end
        ax.errorbar(q_vals, [mean(x) for x in I_coop_freqs],
        [std(x)/sqrt(num_trials) for x in I_coop_freqs],
        label="Q = $Q")
    end
    ax.set_xlim([0,1])
    ax.set_ylim([0,1])
    if ri == 1
        ax.legend(loc=3)
        ax.set_ylabel("coop frequency")
    end
    ax.set_xlabel("q")
    ax.set_title("$rep_norm")
end
plt.tight_layout()

display(fig)
