using Revise
using Institutions
using Statistics
using PyPlot

N = 50
# Q = 2
# q = 0.5

b = 5.0
c = 1.0
w = 1.0
u_s = .025
u_p = .02
u_a = .02

rep_norm = "scoring"

Q_vals = [1, 2, 50]
q_vals = [0.1, 0.3, 0.5, 0.7, 0.9]
coop_means = zeros(Float64, length(Q_vals), length(q_vals))
coop_errors = zeros(Float64, length(Q_vals), length(q_vals))

burn_time = 10000



trial_dist = 10000
num_trials = 1

for (Qi, Q) in enumerate(Q_vals)
    for (qi, q) in enumerate(q_vals)
        #println("simulating Q = $Q, q = $q")
        coop_rate = []
        for i in 1:num_trials
            game = Game(b, c, w, u_s, u_p, u_a)
            pop = Population(N, Q, q, game, rep_norm, [1,2,3,4])
            for j in 1:trial_dist
                evolve!(pop)
            end
            push!(coop_rate, 1.0*sum(pop.prev_actions)/(N*(N-1)))
            println("$Q $q frequencies = $(get_freqs(pop))")
            println("$Q $q reputations = $(get_strat_reputations(pop))")
        end
        coop_means[Qi,qi] = mean(coop_rate)
        coop_errors[Qi,qi] = 2*std(coop_rate)/sqrt(num_trials)
    end
end

fig = plt.figure()
for (Qi, Q) in enumerate(Q_vals)
    plt.errorbar(q_vals, coop_means[Qi,:], coop_errors[Qi,:], label = "Q = $Q")
end
plt.xlim([0,1])
plt.xlabel("q")
plt.ylim([0,1])
plt.ylabel("cooperation frequency")
plt.legend(loc=4)
plt.tight_layout()
display(fig)
