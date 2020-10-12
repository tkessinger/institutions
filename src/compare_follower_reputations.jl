using Revise
using Institutions
using Statistics
using PyPlot

N = 50

num_trials = 1

b = 5.0
c = 1.0

#u_s = 0.0
u_f = 0.02
u_s = 0.025
u_p = 0.05
u_a = 0.05

ϵ = (1 - u_p)*(1 - u_a) + u_p*u_a

strategies = [3]

verbose = false

w = 0.0

inst_type = "strict"

inst_type == "strict" ? q = 0.1 : q = 0.9
Q = 2

E = 1.0
num_gens = 10000
#rep_norm = "simple standing"
rep_norms = ["stern judging", "simple standing", "scoring", "shunning"]
rep_norm = rep_norms[1]

game = Game(b, c, w, u_s, u_p, u_a, "pc")

successful = false
trials = 0
fig, axs = plt.subplots(1,length(rep_norms), figsize = (20, 5), sharey = "row")
for (ri, rep_norm) in enumerate(rep_norms)

        follower_reputations = []
        indep_reputations = []
        follower_freqs = []
        G_vals = []
        g_vals = []
        gn_vals = []
        gf_vals = []

        if rep_norm ∈ ["stern judging", "simple standing"]
                theoretical_rep_value = (1 - u_a)/(2 - ϵ - u_a)
        else
                theoretical_rep_value = u_a/(1 - ϵ + u_a)
        end
        pop = fixation_population(N, E, Q, q, game, rep_norm, strategies, verbose)
        pop.priv_reputations = ones(pop.N,pop.N)
        pop.is_follower[rand(1:N,N÷2)] .= 1
        for i in 1:num_gens
                if rand() < u_f
                        x = rand(1:N)
                        pop.is_follower[x] = 1 - pop.is_follower[x]
                end
                evolve!(pop)
                g = sum(pop.priv_reputations)/pop.N^2
                if rep_norm ∈ ["stern judging", "simple standing"]
                        gn = g*ϵ + (1 - g)*(1 - u_a)
                else
                        gn = g*ϵ + (1 - g)*u_a
                end
                if inst_type == "strict"
                        gf = g^2*(g*ϵ+ (1-g)*(1 - ϵ)) + (1-g^2)*(g*u_a+ (1-g)*(1-u_a))
                else
                        gf = (2*g - g^2)*(g*ϵ+ (1-g)*(1 - ϵ)) + (1 - 2*g + g^2)*(g*u_a+ (1-g)*(1-u_a))
                end
                push!(g_vals, g)
                push!(gf_vals, gf)
                push!(gn_vals, gn)
                push!(G_vals, sum(pop.pub_reputations)/pop.N)
                push!(follower_freqs, 1.0*sum(pop.is_follower)/pop.N)
                push!(follower_reputations, sum(pop.priv_reputations[:, pop.is_follower .== 1])/(pop.N*sum(pop.is_follower)))
                push!(indep_reputations, sum(pop.priv_reputations[:, pop.is_follower .== 0])/(pop.N*sum(1.0 .- pop.is_follower)))
        end
        ax = axs[ri]
        ax.plot(indep_reputations, label = "independent reputation")
        ax.plot(follower_reputations, label = "follower reputation")
        ax.plot(gn_vals, label = "predicted independent reputation")
        ax.hlines(theoretical_rep_value, 0, num_gens, linestyle = "--", label="theory")
        ax.set_title("$rep_norm")
        ax.set_ylim([0, 1])
        if ri == 1
                ax.legend(loc = 4)
        end
end
plt.tight_layout()
display(fig)
# plt.plot(indep_reputations, label = "independent reputation")
# #plt.plot(follower_freqs, ls = "--", label = "follower frequency")
# #plt.plot(gf_vals, label = "predicted follower reputation")
# plt.plot(gn_vals, label = "predicted independent reputation")
# plt.hlines(theoretical_rep_value, 0, num_gens, linestyle = "--")
# plt.ylim([0,1])
# plt.legend(loc = 4)
# plt.tight_layout()
# display(fig)
