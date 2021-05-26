#!/usr/bin/env julia

## parallel_fixation_variable_b.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Parallelized implementation of Institutions.
## Records fixation for institutional adherence
## starting from a single adherent.
## This output corresponds to supplementary figure 14.

using Random, Statistics
using Institutions
using Distributed
using Revise
using ArgParse
using CSV
using Dates
using DataFrames
import JSON

function read_parameters(defpars::Dict{String,Any},
    inputfile = nothing)
	# read and parse JSON file
	# to pass parameters to a worker

    pars = copy(defpars)

    # read JSON file
    if inputfile != nothing
        inpars = JSON.parsefile(inputfile)
    else
        inpars = Dict()
    end

    for parkey in keys(defpars)
        if "type" in keys(pars[parkey])
            if isprimitivetype(pars[parkey]["type"]) ||
                pars[parkey]["type"] == String
                T = pars[parkey]["type"]
            end
        else
            # default type is Float64
            T = Float64
        end

        if T <: Int
            convertf = (val)->round(T, val)
        else
            convertf = (val)->convert(T, val)
        end

        # use defpars for list of usable parameters in JSON
        if parkey in keys(inpars)
            if "value" in keys(inpars[parkey])
                val = inpars[parkey]["value"]
            elseif "range" in keys(inpars[parkey])
                valr = inpars[parkey]["range"]
                if "log" in keys(valr)
                    b = valr["log"]
                    rf = (r)->b.^r
                    pop!(valr, "log")
                else
                    rf = (r)->r
                end
                start = pop!(valr, "start")
                rkws = Dict(zip(Symbol.(keys(valr)), values(valr)))
                val = rf(range(start; rkws...))
            end
        else
            val = pars[parkey]["value"]
        end

        if !isstructtype(typeof(val)) || typeof(val) == String || typeof(val) == Bool
            pars[parkey] = [convertf(val)]
        else
            pars[parkey] = convertf.(val)
        end
    end

    return pars
end

function main(args)
	# main simulation function:
	# parse parameters, take their Cartesian product,
	# run the actual simulation,
	# and record output

    s = ArgParseSettings(description =
        "run ReputationSets simulations across multiple cores")
    @add_arg_table s begin
        "--ncpus"
            arg_type = Int64
            default = max(round(Int, Sys.CPU_THREADS), 1)
        "--input"
            default = nothing
        #"--output"
        #    default=nothing
    end
    parsed_args = parse_args(args, s)

    defpars = Dict{String,Any}([
        "N"     => Dict("value" => 50, "type" => Int64),
        "E"     => Dict("value" => 1.0, "type" => Float64),
		"Q"		=> Dict("value" => 1, "type" => Int64),
		"q"		=> Dict("value" => 1, "type" => Float64),
		"b"     => Dict("value" => 1.0,     "type" => Float64),
		"c"     => Dict("value" => 0.1,     "type" => Float64),
		"w"     => Dict("value" => 1.0,     "type" => Float64),
		"u_s"     => Dict("value" => 0.0,     "type" => Float64),
		"u_p"     => Dict("value" => 0.01,     "type" => Float64),
		"u_a"     => Dict("value" => 0.01,     "type" => Float64),
		"run_type" => Dict("value" => "middle", "type" => String),
		"reputation_norm" => Dict("value" => "stern judging", "type" => String),
        "permitted_strategies" => Dict("value" => "all", "type" => String),
		"independent_board"     => Dict("value" => true,     "type" => Bool),
        "num_trials" => Dict("value" => 50, "type" => Int64),
		"runs_per_trial" => Dict("value" => 50, "type" => Int64),
        "output" => Dict("value" => "output/test.csv", "type" => String)
    ])
    pars = read_parameters(defpars, parsed_args["input"])

    # take the Cartesian product of all parameter combinations
    parsets = collect(Base.product(values(pars)...))
    nsets = length(parsets)

    # setup workers assuming directory is manually added to LOAD_PATH
    addprocs(min(parsed_args["ncpus"], round(Int64, Sys.CPU_THREADS)))
    wpool = WorkerPool(workers())
    #extradir = filter((p)->match(r"/", p) !== nothing, LOAD_PATH)[1]
    extradir = filter((p)->match(r"/", p) !== nothing, LOAD_PATH)
    #@everywhere workers() push!(LOAD_PATH, $extradir)
    [@everywhere workers() push!(LOAD_PATH, $x) for x in extradir]
	# make relevant modules available to each worker
    @everywhere workers() eval(:(using Random))
    @everywhere workers() eval(:(using Statistics))
    @everywhere workers() eval(:(using Institutions))
    @everywhere workers() eval(:(using Dates))

    inputs  = RemoteChannel(()->Channel{Dict}(2 * nsets * maximum(pars["num_trials"])))
    results = RemoteChannel(()->Channel{Dict}(2 * nsets * maximum(pars["num_trials"])))

    @everywhere function run_worker(inputs, results)
        # save trial number and random seed
        seed = Dict(zip(["seed1", "seed2", "seed3", "seed4"], Random.MersenneTwister().seed))

        while true
			# append seed to parameter dictionary
            pard = take!(inputs)
            pard = merge(pard, seed)

            N = pard["N"] # population size
			E = pard["E"] # empathy
			Q = pard["Q"] # institution size
			q = pard["q"] # institution threshold

			reputation_norm = pard["reputation_norm"] # social norm

			# parse permitted strategies
			if pard["permitted_strategies"] == "all"
				permitted_strategies = [1,2,3,4]
			else
				permitted_strategies = parse.(Int64, split(pard["permitted_strategies"], ""))
			end

			# game parameters
			b = pard["b"] # benefit to cooperation
			c = pard["c"] # cost to cooperation
			w = pard["w"] # selection strength
			# is the board external or not?
			independent_board = pard["independent_board"]

			u_s = pard["u_s"] # mutation rate
			u_p = pard["u_p"] # performance error (e1)
			u_a = pard["u_a"] # assessment error (e2)

			# type of run
			run_type = pard["run_type"]

			# num_trials independent runs will ultimately be initialized
			# each one will sample runs_per_trial times
			num_trials = pard["num_trials"]
			runs_per_trial = pard["runs_per_trial"]

			# name of output file
            output = pard["output"]

            println("--- running ", pard["nrun"], " --- ")
            flush(stdout)

			# initialize game
			game = Game(b, c, w, u_s, u_p, u_a, "pc")

			total_interactions = N^2

			# initialize success/failure counters and reputations
			successes = 0
			failures = 0
			equilibrium_reputations = Float64[]
			equilibrium_inst_reputations = Float64[]

			for g in 1:runs_per_trial
				# initialize population with all Disc
				pop = fixation_population(N, E, Q, q, game, reputation_norm, permitted_strategies, independent_board)
				# if multiple strategies are allowed,
				# decide where on the simplex the population begins
				# i.e., alter strategies
				if run_type == "equal"
					pop.strategies[1:(pop.N÷3)] .= 1
					pop.strategies[(pop.N÷3 + 1):(2*pop.N÷3)] .= 2
					pop.strategies[(2* pop.N÷3 + 1):pop.N] .= 3
				elseif run_type == "AllC"
					pop.strategies[1:(pop.N-4)] .= 1
					pop.strategies[(pop.N-3):(pop.N-2)] .= 2
					pop.strategies[(pop.N-1):pop.N] .= 3
				elseif run_type == "AllD"
					pop.strategies[1:2] .= 1
					pop.strategies[3:(pop.N-2)] .= 2
					pop.strategies[(pop.N-1):pop.N] .= 3
				elseif run_type == "Disc"
					pop.strategies[1:2] .= 1
					pop.strategies[3:4] .= 2
					pop.strategies[5:pop.N] .= 3
				else
				end
				# allow reputations to equilibrate
				for i in 1:100
					update_actions_and_fitnesses!(pop)
					update_reputations!(pop)
				end
				# record initial reputations
				init_reputation = sum(pop.priv_reputations)/N^2
				init_inst_reputation = sum(pop.pub_reputations)/N
				push!(equilibrium_reputations, init_reputation)
				push!(equilibrium_inst_reputations, init_inst_reputation)
				# if a discriminator exists, randomly make one a non-follower
				# else, make a discriminator and make it a non-follower
				pop.is_follower .= 1
				if sum([pop.strategies[x] .== 3 for x in 1:pop.N]) > 0
					pop.is_follower[rand(filter(x -> pop.strategies[x] .== 3, 1:pop.N))] = 0
				else
					x = rand(1:pop.N)
					pop.strategies[x] = 3
					pop.is_follower[x] = 0
				end
				# evolve population until fixation or extinction
				while sum(pop.is_follower) ∉ [0, pop.N]
					evolve!(pop)
				end
				# record success (fixation) or failure (extinction)
				if sum(pop.is_follower) == pop.N
					successes += 1
				else
					failures += 1
				end
			end

			# record the average reputation over many runs
			reputations = mean(equilibrium_reputations)
			inst_reputations = mean(equilibrium_inst_reputations)

			pard["successes"] = successes
			pard["failures"] = failures
			pard["reputations"] = reputations
			pard["inst_reputations"] = inst_reputations

            # return data to master process
            put!(results, pard)
        end
    end

    total_time_start = now()

    # load parameter sets into inputs channel
    nruns = 0
    for parset in parsets
        pard = Dict(zip(keys(pars), parset))
        println("--- queueing --- ")
        foreach(k->print(k, ": ", pard[k], ", "), sort(collect(keys(pard))))
        println()
        flush(stdout)
        for rep in 1:pard["num_trials"]
            nruns += 1
            rpard = copy(pard)
            rpard["rep"] = rep
            rpard["nrun"] = nruns
            put!(inputs, rpard)
        end
    end

    # start workers running on parameter sets in inputs
    for w in workers() # start tasks on the workers to process requests in parallel
        remote_do(run_worker, w, inputs, results)
    end

    # create output file name and data table
    output = pars["output"][1]
    println(output)
    file = occursin(r"\.csv$", output) ? output : output * ".csv"
    cols = push!(sort(collect(keys(pars))),
                 ["rep", "successes", "failures", "reputations", "inst_reputations", "seed1", "seed2", "seed3", "seed4"]...)
    dat = DataFrame(Dict([(c, Any[]) for c in cols]))

    # grab results and output to CSV
    for sim in 1:nruns
        # get results from parallel jobs
        flush(stdout)
        resd = take!(results)
        nrun = pop!(resd, "nrun")

        # add to table (must convert dict keys to symbols) and save
        push!(dat, Dict([(Symbol(k), resd[k]) for k in keys(resd)]))
        CSV.write(file, dat)
    end
    total_time_stop = now()
    total_time = Dates.canonicalize(Dates.CompoundPeriod(round(total_time_stop - total_time_start, Dates.Second(1))))
    println("total time elapsed: $total_time")
end

# specify input file here
main(["--input", "submit/fixation_paper_adherents.json"])
