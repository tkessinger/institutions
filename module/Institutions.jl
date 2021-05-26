#!/usr/bin/env julia

## Institutions.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Module for simulating cooperation and defection
## when reputations are tracked by institutions.
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

module Institutions

	using Random, StatsBase, Combinatorics, Distributions

	export Game, Population
	export empathy_population, institution_population, fixation_population
	export evolve!
	export get_freqs, get_pub_reputations, get_priv_reputations
	export get_strat_pub_reputations, get_strat_priv_reputations
	# following exports are for test purposes
	export update_strategies_pc!
	export update_reputations!, update_actions_and_fitnesses!
	export reputation_norm, determine_action, mutate!

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
		Q::Int64 # not mutable: institution size
		q::Float64 # not mutable: threshold of institution members
			# who need to agree that someone's reputation is ``good''
		game::Game # Game object, see struct Game above
		norm::String # reputation norm
		strategies::Array{Int64, 1} # array of reputation assessment strategies
		empathies::Array{Float64, 1} # individual empathy values
		is_follower::Array{Bool, 1} # does this person follow an institution or not?
		priv_reputations::Array{Int64, 2} # private assessments
			# of the entire population
		pub_reputations::Array{Int64, 1} # institution's public reputation
		prev_actions::Array{Int64, 2} # the last action each individual took toward each other
		fitnesses::Array{Float64, 1} # array of fitnesses
		permitted_strategies::Array{Int64, 1} # which strategies are allowed to appear
		generation::Int64 # current generation
		independent_board::Bool # make the institution independent or not?
		verbose::Bool # turn this on for error tracking

	end

	function empathy_population(
		N::Int64,
		E::Float64,
		game::Game,
		norm::String,
		initial_strategies::Array{Int64, 1}=[1,2,3,4],
		independent_board::Bool=false
		)
		# for generating a population of empathetic individuals
		# who use private assessment only
		# see Radzvilavicius and Plotkin (2019) for empathy explanation
		# begin by initializing the population with random strategies
		strategies = rand(initial_strategies, N)
		Q = 0 # no institution
		q = 0.0
		empathies = E*ones(Float64, N)
		is_follower = zeros(Bool, N)
		priv_reputations = rand([0,1], N, N)
		pub_reputations = zeros(Int64, N)
		prev_actions = zeros(Int64, N, N)
		fitnesses = zeros(Float64, N)
		permitted_strategies = initial_strategies
		generation = 0
		return Population(N, Q, q, game, norm, strategies, empathies,
			is_follower, priv_reputations, pub_reputations,
			prev_actions, fitnesses, permitted_strategies,
			generation, independent_board, false)
	end

	function institution_population(
		N::Int64,
		Q::Int64,
		q::Float64,
		game::Game,
		norm::String,
		initial_strategies::Array{Int64, 1}=[1,2,3,4],
		independent_board::Bool=false
		)
		# for generating a population of institution adherents
		# begin by initializing the population with random strategies
		strategies = rand(initial_strategies, N)
		empathies = zeros(Float64, N)
		is_follower = ones(Bool, N)
		priv_reputations = rand([0,1], N, N)
		pub_reputations = zeros(Int64, N)
		[pub_reputations[x] = (sum(priv_reputations[1:Q,x]) > q*Q) for x in 1:N]
		prev_actions = zeros(Int64, N, N)
		fitnesses = zeros(Float64, N)
		permitted_strategies = initial_strategies
		generation = 0
		return Population(N, Q, q, game, norm, strategies, empathies,
			is_follower, priv_reputations, pub_reputations,
			prev_actions, fitnesses, permitted_strategies,
			generation, independent_board, false)
	end

	function fixation_population(
		N::Int64,
		E::Float64,
		Q::Int64,
		q::Float64,
		game::Game,
		norm::String,
		initial_strategies::Array{Int64, 1}=[3],
		independent_board::Bool=false
		)
		# for studying the invasion of institution adherence
		# in a population of (potentially empathetic) non-adherents
		# begin by initializing the population with random strategies
		strategies = rand(initial_strategies, N)
		empathies = E*ones(Float64, N)
		is_follower = zeros(Bool, N)
		priv_reputations = rand([0,1], N, N)
		pub_reputations = zeros(Int64, N)
		[pub_reputations[x] = (sum(priv_reputations[1:Q,x]) > q*Q) for x in 1:N]
		prev_actions = zeros(Int64, N, N)
		fitnesses = zeros(Float64, N)
		permitted_strategies = initial_strategies
		generation = 0
		return Population(N, Q, q, game, norm, strategies, empathies,
			is_follower, priv_reputations, pub_reputations,
			prev_actions, fitnesses, permitted_strategies,
			generation, independent_board, false)
	end

	function get_freqs(
		pop::Population
		)
		# return strategy frequencies
		return [1.0*sum(pop.strategies .== x)/pop.N for x in 1:4]
	end

	function get_pub_reputations(
		pop::Population
		)
		# return fraction G of good public reputations
		return 1.0*sum(pop.pub_reputations .== 1)/pop.N
	end

	function get_priv_reputations(
		pop::Population
		)
		# return fraction g of good private reputations
		return 1.0*sum(pop.priv_reputations .== 1)/pop.N^2
	end

	function get_strat_pub_reputations(
		pop::Population
		)
		# return G_i (by strategy type)
		strat_reps = [1.0*sum((pop.strategies .== x) .* (pop.pub_reputations .== 1))/max(1,sum(pop.strategies .== x)) for x in 1:4]
		return strat_reps
	end

	function get_strat_priv_reputations(
		pop::Population
		)
		# return g_i (by strategy type)
		strat_reps = [sum(sum([(pop.strategies .== x) .* (pop.priv_reputations[n,:] .== 1)
			for n in 1:pop.N]))/max(1,sum(pop.N*(pop.strategies .== x))) for x in 1:4]
		return strat_reps
	end

	function mutate!(
		pop::Population
		)
		# randomly mutate a single individual to a new strategy
		if rand() < pop.game.u_s
			n = rand(1:pop.N)
			# this line allows the following mutations:
			# AllC or AllD to Disc or RDisc
			# Disc or RDisc to AllC or AllD
			# this comes from a model where strategies are represented
			# as [p, q]
			# with p the probability of cooperating with good
			# and q the probability of cooperating with bad
			# AllC = [1,1], AllD = [0,0], Disc = [1,0], RDisc = [0,1]
			# and only one "step" is allowed.
			# changing this to allow "free" mutation has minimal effect.
			old_strat = pop.strategies[n]
			new_strat = rand(pop.permitted_strategies)
			#new_strat = ((4 - pop.strategies[n])÷2)*2+rand([1,2])
			# if the new strategy is permitted, use it
			if new_strat ∈ pop.permitted_strategies
				pop.strategies[n] = new_strat
				if pop.verbose println("mutated $n from $old_strat to $new_strat") end
			end
		end
	end

	function update_strategies_pc!(
		pop::Population
		)
		# chooses a random pair of individuals to compare via a sigmoid function
		# the fitter individual has a chance of invading the less fit one
		# (i.e., forcing them to change strategy)
		# invader, invadee = sample(1:pop.N, 2)
		invader, invadee = rand(1:pop.N, 2) # this allows the same individual to be picked for both
		# sigmoid update function
		# sanity check: this should be higher if invader fitness > invadee fitness
		update_function = 1.0/(1.0+exp(-pop.game.w*(pop.fitnesses[invader]-pop.fitnesses[invadee])))
		if pop.verbose println("randomly chosen invader $invader and invadee $invadee") end
		if pop.verbose println("fitnesses are $(pop.fitnesses[invader]) and $(pop.fitnesses[invadee])") end
		if pop.verbose println("strategies are $(pop.strategies[invader]) and $(pop.strategies[invadee])") end
		if pop.verbose println("update function is $update_function") end
		# if rand() returns a number smaller than the update value,
		# update strategy, empathy value, and follower status
		if rand() < update_function
			pop.strategies[invadee] = pop.strategies[invader]
			pop.is_follower[invadee] = pop.is_follower[invader]
			pop.empathies[invadee] = pop.empathies[invader]
			if pop.verbose println("$invadee adopts strategy $(pop.strategies[invadee])") end
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
			mutate!(pop)
			update_actions_and_fitnesses!(pop)
			# then make sure everyone's reputations are updated
			if pop.verbose println("updating reputations") end
			update_reputations!(pop)
			# then, finally, select a pair of individuals whose fitnesses we will compare
			if pop.verbose println("updating strategy") end
			# only pairwise comparison is currently supported
			update_strategies_pc!(pop)
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

		# for each donor i
		for i in 1:pop.N
			# for each recipient j
			for j in 1:pop.N
				# i behaves toward j based on i's view of j
				# (note: i may already have adopted the institutional view)
				if pop.is_follower[i]
					i_intended_action = determine_action(pop.strategies[i], pop.pub_reputations[j])
				else
					i_intended_action = determine_action(pop.strategies[i], pop.priv_reputations[i,j])
				end
				i_rand = rand()
				i_rand > pop.game.u_p ? i_action = i_intended_action : i_action = 0
				new_actions[i,j] = i_action
			end
		end
		for i in 1:pop.N
			# i receives benefit b for each time it is cooperated with
			new_fitnesses[i] += pop.game.b*sum(new_actions[:,i])/pop.N
			# i pays cost c for each time it cooperates
			new_fitnesses[i] -= pop.game.c*sum(new_actions[i,:])/pop.N
			# note: these include self-interactions, which cannot convey net cost or benefit
		end

		if pop.verbose println("new fitnesses look like $new_fitnesses") end
		if pop.verbose println("new actions look like $new_actions") end

		# assign the new fitnesses and actions
		pop.fitnesses = new_fitnesses
		pop.prev_actions = new_actions
	end

	function draw_reputation(
		pop::Population,
		i::Int64,
		j::Int64,
		k::Int64,
		board::Bool = false
		)
		# performs a single "draw" of reputation for a given observer, donor, and recipient
		# the "board" parameter is passed in as true IFF we are looking at the
		# institutional assessments of an external board,
		# who presumably do not care about empathy at all.
		# they evaluate reputation based solely on
		# adherence to the institution's assessments.
		# otherwise it is passed in as false.

		# examine what j did to k
		action = pop.prev_actions[j,k]
		# if empathetic...
		if rand() < pop.empathies[i]
			# if j is a follower, look at the institutional view of k
			if pop.is_follower[j] || board == true
				normed_reputation = reputation_norm(action, pop.pub_reputations[k], pop.norm)
			# otherwise, consider j's private view of k
			else
				normed_reputation = reputation_norm(action, pop.priv_reputations[j,k], pop.norm)
			end
		# if not empathetic...
		else
			# if i is a follower, look at the institutional view of k
			if pop.is_follower[i] || board == true
				normed_reputation = reputation_norm(action, pop.pub_reputations[k], pop.norm)
			# otherwise, consider i's private view of k
			else
				normed_reputation = reputation_norm(action, pop.priv_reputations[i,k], pop.norm)
			end
		end
		return normed_reputation
	end

	function update_reputations!(
		pop::Population
		)
		# update each individual's public reputation
		# by choosing a random action to observe
		# then adjust each individual's private attitude about every other individual

		# initialize public and private reputations at zero
		new_priv_reputations = zeros(Int64, pop.N, pop.N)
		new_pub_reputations = zeros(Int64, pop.N)
		# for each observer
		for i in 1:pop.N
			# for each donor
			for j in 1:pop.N
				# check a random other individual k and see what j did to k
				k = rand(1:pop.N)
				normed_reputation = draw_reputation(pop, i, j, k)
				# apply assessment error
				rep_rand = rand()
				rep_rand > pop.game.u_a ? new_priv_reputations[i,j] = normed_reputation : new_priv_reputations[i,j] = 1 - normed_reputation
				# if pop.verbose
				# 	println("$i analyzes $j's behavior toward $k")
				# 	println("$j did $action and $k's reputation is $(pop.pub_reputations[k])")
				# 	println("the norm \"$(pop.norm)\" predicts that $j's reputation be $normed_reputation")
				# 	println("the random number is $rep_rand")
				# 	println("$i judges $j's reputation to be $(new_priv_reputations[i,j])")
				# end
			end
		end
		# if an individual's total reputation, summed over institution members,
		# is greater than q*Q, their reputation is broadcast as good
		# otherwise, it is broadcast as bad
		# initialize the "dummy" set of reputations for institution members
		dummy_reputations = zeros(Int64, pop.Q, pop.N)
		# for each observer (institution member)
		for i in 1:pop.Q
			# for each donor
			for j in 1:pop.N
				# check a random other individual k and see what j did to k
				k = rand(1:pop.N)
				# if the board is external, pop.independent board will override empathy
				normed_reputation = draw_reputation(pop, i, j, k, pop.independent_board)
				rep_rand = rand()
				# apply assessment error
				rep_rand > pop.game.u_a ? dummy_reputations[i,j] = normed_reputation : dummy_reputations[i,j] = 1 - normed_reputation
			end
		end
		for j in 1:pop.N
			new_pub_reputations[j] = (sum(dummy_reputations[:,j]) > pop.q*pop.Q)
		end

		pop.priv_reputations = new_priv_reputations
		pop.pub_reputations = new_pub_reputations
		if pop.verbose
			println("private reputations are $new_priv_reputations and cutoff is $(pop.q*pop.Q)")
			println("reputations are broadcast to be $new_pub_reputations")
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
