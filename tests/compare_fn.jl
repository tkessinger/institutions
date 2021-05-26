using Revise
using Institutions
using Statistics
using PyPlot

N = 500

b = 5.0
c = 1.0

#u_s = 0.0
u_f = 0.0
u_s = 1.0/N
u_p = 1.0/N
u_a = 1.0/N

ϵ = (1 - u_p)*(1 - u_a) + u_p*u_a

strategies = [1,2,3]

verbose = false

w = 1.0

inst_type = "strict"

inst_type == "strict" ? q = 0.9 : q = 0.1
Q = 1

E = 0.0
num_gens = 1000

#rep_norm = "simple standing"
rep_norms = ["stern judging", "simple standing", "scoring", "shunning"]
#rep_norm = rep_norms[1]

game = Game(b, c, w, u_s, u_p, u_a, "pc")


fig, axs = plt.subplots(length(rep_norms), 4, figsize = (20, 5*length(rep_norms)), sharey = "row", sharex = "col")

dg = 10

for (ri, rep_norm) in enumerate(rep_norms)
    strat_freqs = []
    Gi_vals = []
    gi_vals = []
    predicted_Gi_vals = []

    g_vals = []
    G_vals = []

    pop = fixation_population(N, E, Q, q, game, rep_norm, strategies, verbose)
    pop.is_follower[rand(1:N,N÷2)] .= 1


    for i in 1:num_gens
        #update_actions_and_fitnesses!(pop)
        #update_reputations!(pop)
        evolve!(pop)
        g = sum(pop.priv_reputations)/pop.N^2
        g1 = sum(pop.priv_reputations[:, pop.strategies .== 1])/(pop.N*sum(pop.strategies .== 1))
        g2 = sum(pop.priv_reputations[:, pop.strategies .== 2])/(pop.N*sum(pop.strategies .== 2))
        g3 = sum(pop.priv_reputations[:, (pop.strategies .== 3) .& (pop.is_follower .== 0)])/(pop.N*sum((pop.strategies .== 3) .& (pop.is_follower .== 0)))
        g4 = sum(pop.priv_reputations[:, (pop.strategies .== 3) .& (pop.is_follower .== 1)])/(pop.N*sum((pop.strategies .== 3) .& (pop.is_follower .== 1)))

        push!(g_vals, g)
        push!(gi_vals, [g1, g2, g3, g4])

        if inst_type == "strict"
            push!(predicted_Gi_vals, [g1^2, g2^2, g3^2, g4^2])
        else
            push!(predicted_Gi_vals, [1 - (1 - g1)^2, 1 - (1 - g2)^2, 1 - (1 - g3)^2, 1 - (1 - g4)^2])
        end

        G = sum(pop.pub_reputations)/pop.N
        G1 = sum(pop.pub_reputations[pop.strategies .== 1])/(sum(pop.strategies .== 1))
        G2 = sum(pop.pub_reputations[pop.strategies .== 2])/(sum(pop.strategies .== 2))
        G3 = sum(pop.pub_reputations[(pop.strategies .== 3) .& (pop.is_follower .== 0)])/(sum((pop.strategies .== 3) .& (pop.is_follower .== 0)))
        G4 = sum(pop.pub_reputations[(pop.strategies .== 3) .& (pop.is_follower .== 1)])/(sum((pop.strategies .== 3) .& (pop.is_follower .== 1)))

        push!(G_vals, G)
        push!(Gi_vals, [G1, G2, G3, G4])

        push!(strat_freqs, 1.0/pop.N*[sum(pop.strategies .== 1), sum(pop.strategies .== 2), sum((pop.strategies .== 3) .& (pop.is_follower .== 0)), sum((pop.strategies .== 3) .& (pop.is_follower .== 1))])
    end

    strat_freqs = permutedims(hcat(strat_freqs...))
    gi_vals = permutedims(hcat(gi_vals...))
    Gi_vals = permutedims(hcat(Gi_vals...))
    predicted_Gi_vals = permutedims(hcat(predicted_Gi_vals...))

    ax = axs[ri,1]
    [ax.plot(collect(1:num_gens)[1:dg:end], gi_vals[1:dg:end,i], label = "g_$i") for i in 1:4]
    ax.plot(collect(1:num_gens)[1:dg:end], g_vals[1:dg:end], label = "g")
    ax.legend(loc=1)
    ax = axs[ri,2]
    ax.set_title("$rep_norm")
    [ax.plot(collect(1:num_gens)[1:dg:end], Gi_vals[1:dg:end,i], label = "G_$i") for i in 1:4]
    ax.plot(collect(1:num_gens)[1:dg:end], G_vals[1:dg:end], label = "G")
    ax.legend(loc=1)
    ax = axs[ri,3]
    [ax.plot(collect(1:num_gens)[1:dg:end], predicted_Gi_vals[1:dg:end,i], label = "predicted G_$i") for i in 1:4]
    ax.legend(loc=1)
    ax = axs[ri,4]
    [ax.plot(collect(1:num_gens)[1:dg:end], strat_freqs[1:dg:end,i], label = "f_$i") for i in 1:4]
    ax.legend(loc=1)
end

[ax.set_ylim([0,1]) for ax in axs]

plt.tight_layout()

display(fig)
