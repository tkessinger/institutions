#!/usr/bin/env julia

## Institutions.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Module for simulating cooperation and defection
## when reputations are tracked individually
## but individuals can empathetically mimic others' judgments.
## Adapted from work by Arunas Radzvilavicius.

# Outline:
# 1. Initialize population.
# 2. Assign initial private reputations for each institution member.
# 3. Broadcast these reputations to public.
# 4. Compute fitness by having each individual act on each other individual
# based on their reputation and the social norm.
# With small probability u, cooperators accidentally defect.
# 5. Keep track of those actions.
# 6. Update reputations by having each institution member observe a single random (past)
# action of every other individual.
# 7. Broadcast these reputations to public.
# 8. Choose a random individual to update their strategy via a sigmoid function.
# 9. Goto 4 and repeat.

# TODO: add empathy parameter, write new reputations function

module Empathy

	using Random, StatsBase, Combinatorics

	export Game, Population
	export evolve!
	export update_strategies_db!, update_strategies_pc!
	export update_reputations!, update_actions_and_fitnesses!
	export get_freqs, get_reputations, get_strat_reputations
	export reputation_norm

	struct Game
		# structure for storing game parameters and such

		b::Float64 # benefit to cooperating
		c::Float64 # cost to cooperating
		w::Float64 # selection strength
		u_s::Float64 # mutation rate
		u_p::Float64 # probability of choosing the "wrong" action
		u_a::Float64 # probability of incorrectly assigning reputation
		update_rule::String # update rule: PC or death-birth
		A::Array{Float64, 2} # the actual game matrix

		function Game(
			b::Float64,
			c::Float64,
			w::Float64,
			u_s::Float64,
			u_p::Float64,
			u_a::Float64,
			update_rule::String
			)
			return new(b, c, w, u_s, u_p, u_a,
				update_rule, [0.0 b; -c b-c])
		end

		function Game(
			b::Float64,
			c::Float64,
			w::Float64,
			u_s::Float64,
			u_p::Float64,
			u_a::Float64,
			)
			return new(b, c, w, u_s, u_p, u_a,
				"pc", [0.0 b; -c b-c])
		end
	end

	mutable struct Population
		# the population of individuals and all information about them

		N::Int64 # not mutable: population size
		E::Float64 # not mutable: empathy parameter
		game::Game
		norm::String
		strategies::Array{Int64, 1} # array of reputation assessment strategies
		priv_reputations::Array{Int64, 2} # institution members' private assessments
			# of the entire population
		prev_actions::Array{Int64, 2} # the last action each individual took toward each other
		fitnesses::Array{Float64, 1} # array of fitnesses
		permitted_strategies::Array{Int64, 1} # which strategies are allowed to appear
		generation::Int64 # current generation
		verbose::Bool # turn this on for error tracking

		# constructor if sets and game are already specified
		function Population(
			N::Int64,
			E::Float64,
			game::Game,
			norm::String,
			initial_strategies::Array{Int64, 1}=[1,2,3,4],
			verbose::Bool=false
			)
			# begin by initializing the population with random strategies
			strategies = rand(initial_strategies, N)
			#strategies = 2*ones(Int64, N)
			# randomize individuals' private reputations
			priv_reputations = rand([0, 1], N, N)
			# previous actions and fitnesses start at zero (these will get updated)
			prev_actions = zeros(Int64, N, N)
			fitnesses = zeros(Float64, N)
			permitted_strategies = initial_strategies
			generation = 0
			return new(N, E, game, norm, strategies, priv_reputations,
				prev_actions, fitnesses, permitted_strategies, generation, verbose)
		end
	end

	function get_freqs(
		pop::Population
		)
		return [1.0*sum(pop.strategies .== x)/pop.N for x in 1:4]
	end

	function get_reputations(
		pop::Population
		)
		return 1.0*sum(pop.pub_reputations .== 1)/pop.N
	end

	function get_strat_reputations(
		pop::Population
		)
		strat_reps = [1.0*sum((pop.strategies .== x) .* (pop.pub_reputations .== 1))/max(1,sum(pop.strategies .== x)) for x in 1:4]
		return strat_reps
	end


	function update_strategies_db!(
		pop::Population
		)
		# randomly choose someone to die
		invadee = sample(1:pop.sets.N)
		# compute the fitnesses of every other individual in the population
		invasion_fitnesses = 1.0 .- pop.game.w .+ pop.game.w*pop.fitnesses[filter(x->x!=invadee, 1:pop.sets.N)]
		# choose a random other individual in the population, weighted by fitness0
		invader = sample(filter(x->x!=invadee, collect(1:pop.sets.N)), Weights(invasion_fitnesses))
		# the chosen individual's strategy replaces the deceased
		if pop.verbose println("randomly chosen invader $invader and invadee $invadee") end
		if pop.verbose println("fitnesses are $(pop.fitnesses[invader]) and $(pop.fitnesses[invadee])") end
		if pop.verbose println("all fitnesses are $(pop.fitnesses)") end
		if rand() < pop.game.u_s
			# this is where we allow the invadee to mutate
			#pop.strategies[invadee] = rand(filter(x->x!=pop.strategies[invader], pop.permitted_strategies))
			pop.strategies[invadee] = rand(pop.permitted_strategies)
			if pop.verbose println("mutating $invadee to strategy $(pop.strategies[invadee])") end
		else
			pop.strategies[invadee] = pop.strategies[invader]
			if pop.verbose println("adopting strategy $(pop.strategies[invadee])") end
		end
	end

	function update_strategies_pc!(
		pop::Population
		)
		# chooses a random pair of individuals to compare via a sigmoid function
		# the fitter individual has a chance of invading the less fit one
		# (i.e., forcing them to change strategy)
		invader, invadee = sample(1:pop.N, 2)
		# sigmoid update function
		# sanity check: this should be higher if invader fitness > invadee fitness
		update_function = 1.0/(1.0+exp(-pop.game.w*(pop.fitnesses[invader]-pop.fitnesses[invadee])))
		if pop.verbose println("randomly chosen invader $invader and invadee $invadee") end
		if pop.verbose println("fitnesses are $(pop.fitnesses[invader]) and $(pop.fitnesses[invadee])") end
		if pop.verbose println("strategies are $(pop.strategies[invader]) and $(pop.strategies[invadee])") end
		if pop.verbose println("update function is $update_function") end
		if rand() < update_function
			pop.strategies[invadee] = pop.strategies[invader]
			if pop.verbose println("$invadee adopts strategy $(pop.strategies[invadee])") end
		end
		if rand() < pop.game.u_s
			# this is where we allow the invadee to mutate
			pop.strategies[invadee] = rand(pop.permitted_strategies)
			if pop.verbose println("mutating $invadee to strategy $(pop.strategies[invadee])") end
		end
	end

	function evolve!(
		pop::Population,
		generations::Int64 = 1
		)
		# the main evolution function
		# generations allows us to specify how many generations to let the simulation run for
		for i in 1:generations
			if pop.verbose println("initiating generation $(pop.generation)") end
			# we first need to choose actions and update fitnesses
			if pop.verbose println("updating actions and fitnesses") end
			update_actions_and_fitnesses!(pop)
			# then make sure everyone's reputations are updated
			if pop.verbose println("updating reputations") end
			update_reputations_empathy!(pop)
			# then, finally, select a pair of individuals whose fitnesses we will compare
			if pop.verbose println("updating strategy") end
			if pop.game.update_rule ∈ ["pc", "pairwise_comparison", "im", "imitation"]
				update_strategies_pc!(pop)
			elseif pop.game.update_rule ∈ ["db", "death_birth"]
				update_strategies_db!(pop)
			end
			pop.generation += 1
		end
	end

	function determine_action(
		strategy::Int64,
		recipient_reputation::Int64
		)
		# for a strategy and reputation, spits out an action
		if strategy == 1
			# if AllC, always cooperate
			return 1
		elseif strategy == 2
			# if AllD, always defect
			return 0
		elseif strategy == 3
			# if Disc, behave according to target's reputation
			return recipient_reputation
		elseif strategy == 4
			# if RDisc, behave according to opposite of target's reputation
			return 1 - recipient_reputation
		end
	end

	function update_actions_and_fitnesses!(
		pop::Population
		)
		# choose actions between every pair of individuals in every set,
		# then update their fitnesses accordingly

		# initialize all fitnesses and actions at zero
		new_fitnesses = zeros(Float64, pop.N)
		new_actions = zeros(Int64, pop.N, pop.N)
		# for each pair of individuals
		for (i, j) in filter(x -> x[1] < x[2], collect(Base.product(1:pop.N,1:pop.N)))
			if pop.verbose println("updating actions of $i and $j") end

			# if the random number is larger than the error rate, do the intended action
			# else, defect
			i_rand = rand()
			i_intended_action = determine_action(pop.strategies[i], pop.priv_reputations[i,j])
			i_rand > pop.game.u_p ? i_action = i_intended_action : i_action = 0
			# repeat with j
			j_rand = rand()
			j_intended_action = determine_action(pop.strategies[j], pop.priv_reputations[j,i])
			j_rand > pop.game.u_p ? j_action = j_intended_action : j_action = 0

			if pop.verbose
				println("$i's strategy is $(pop.strategies[i])")
				println("$j's reputation is $(pop.priv_reputations[j,i])")
				println("$i's random number was $i_rand")
				println("$i intended to do $i_intended_action and did $i_action")
				println("$j's strategy is $(pop.strategies[j])")
				println("$i's reputation is $(pop.priv_reputations[i,j])")
				println("$j's random number was $j_rand")
				println("$j intended to do $j_intended_action and did $j_action")
				println("$i earns a payoff of $(pop.game.A[i_action+1, j_action+1])")
				println("$j earns a payoff of $(pop.game.A[j_action+1, i_action+1])")
			end
			# adjust i and j's fitnesses according to their strategies and the game matrix
			new_fitnesses[i] += pop.game.A[i_action+1, j_action+1] # +1 because julia is 1-indexed
			new_fitnesses[j] += pop.game.A[j_action+1, i_action+1]
			# store their last actions toward each other in this set
			new_actions[i,j] = i_action
			new_actions[j,i] = j_action
		end
		if pop.verbose println("new fitnesses look like $new_fitnesses") end
		if pop.verbose println("new actions look like $new_actions") end

		pop.fitnesses = new_fitnesses
		pop.prev_actions = new_actions
	end

	function update_reputations_empathy!(
		pop::Population
		)
		# update each individual's public reputation within each set
		# by choosing a random action to observe
		# then adjust each individual's private attitude about every other individual

		### CONVENTION:
		### reputations[i,j] is
		### i's view of j's reputation!

		# initialize public and private reputations at zero
		new_priv_reputations = zeros(Int64, pop.N, pop.N)
		# convention:
		# i: observer
		# j: donor
		# k: recipient
		# for each individual
		for i in 1:pop.N
			# evaluate the reputation of each donor
			for j in filter(x -> x != i, 1:pop.N)
				# look at what the donor did to a random recipient
				k = rand(filter(x -> x ∉ [i,j], 1:pop.N))
				action = pop.prev_actions[j,k]
				empathy_rand = rand()
				# if empathy_rand is less than E, i uses the donor's (j) view of k
				# else, i uses their own (i) view of k
				empathy_rand < pop.E ? k_reputation = pop.priv_reputations[j,k] : k_reputation = pop.priv_reputations[i,k]
				# apply the norm to determine i's attitude toward j
				# (i.e., j's reputation in i's eyes)
				normed_reputation = reputation_norm(action, k_reputation, pop.norm)
				rep_rand = rand()
				rep_rand > pop.game.u_a ? new_priv_reputations[i,j] = normed_reputation : new_priv_reputations[i,j] = 1 - normed_reputation
				if pop.verbose
					println("$i analyzes $j's behavior toward $k")
					println("$j did $action to $k")
					empathy_rand < pop.E ? println("$i is empathetic and views $k's reputation as $(pop.priv_reputations[j,k])") : println("$i is not empathetic and view's $k's reputation as $(pop.priv_reputations[i,k])")
					println("the norm \"$(pop.norm)\" predicts that $j's reputation be $normed_reputation")
					println("the random number is $rep_rand")
					println("$i judges $j's reputation to be $(new_priv_reputations[i,j])")
				end
			end
		end
		pop.priv_reputations = new_priv_reputations
		# if an individual's total reputation, summed over institution members,
		# is greater than q*Q, their reputation is broadcast as good
		# otherwise, it is broadcast as bad
		if pop.verbose
			println("private reputations are $new_priv_reputations")
		end
	end

	function reputation_norm(
		action::Int64, # action taken toward the recipient
		reputation::Int64, # reputation of recipient
		norm_ID::String # identity of norm used
		)
		if norm_ID == "stern judging"
			# stern judging norm
			# cooperating with G and defecting with B are good
			# anything else is bad
			norm = [1 0; 0 1]
		elseif norm_ID == "shunning"
			# shunning norm
			# cooperating with G is good
			# anything else is bad
			norm = [0 0; 0 1]
		elseif norm_ID == "scoring"
			# scoring norm
			# cooperating is good
			# defecting is bad
			norm = [0 0; 1 1]
		elseif norm_ID == "simple standing"
			# simple standing norm
			# defecting with G is bad
			# anything else is good
			norm = [1 0; 1 1]
		end
		return norm[action+1, reputation+1] # +1 because julia is 1-indexed
	end

end
# final end statement to close the module
