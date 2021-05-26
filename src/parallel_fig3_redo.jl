#!/usr/bin/env julia

## parallel_fig3_redo.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Parallelized implementation of Institutions.
## Initializes a population of institutional adherents,
## allows them to evolve and update strategy,
## and records frequency, reputation, and cooperation level.

using Random, Statistics
using Institutions
using Distributed
using Revise
using ArgParse
using CSV
using Dates
using DataFrames
import JSON

function read_parameters(
	defpars::Dict{String,Any},
    inputfile = nothing
	)
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
        #println(parkey, T)
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
        "Q"     => Dict("value" => 1, "type" => Int64),
		"q"     => Dict("value" => 0.5,     "type" => Float64),
		"b"     => Dict("value" => 1.0,     "type" => Float64),
		"c"     => Dict("value" => 0.1,     "type" => Float64),
		"w"     => Dict("value" => 1.0,     "type" => Float64),
		"u_s"     => Dict("value" => 0.01,     "type" => Float64),
		"u_p"     => Dict("value" => 0.01,     "type" => Float64),
		"u_a"     => Dict("value" => 0.01,     "type" => Float64),
		"reputation_norm" => Dict("value" => "stern judging", "type" => String),
        "permitted_strategies" => Dict("value" => "all", "type" => String),
        "num_gens" => Dict("value" => 10000, "type" => Int64),
        "num_trials" => Dict("value" => 50, "type" => Int64),
        "output" => Dict("value" => "output/test.csv", "type" => String)
    ])
    pars = read_parameters(defpars, parsed_args["input"])

    # take the Cartesian product of all parameter combinations
    parsets = collect(Base.product(values(pars)...))
    nsets = length(parsets)

    # setup workers assuming directory is manually added to LOAD_PATH
    addprocs(min(parsed_args["ncpus"], round(Int64, Sys.CPU_THREADS / 2)))
    wpool = WorkerPool(workers())
    #extradir = filter((p)->match(r"/", p) !== nothing, LOAD_PATH)[1]
    extradir = filter((p)->match(r"/", p) !== nothing, LOAD_PATH)
    #@everywhere workers() push!(LOAD_PATH, $extradir)
    [@everywhere workers() push!(LOAD_PATH, $x) for x in extradir]
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
            pard = take!(inputs)
            pard = merge(pard, seed)

            N = pard["N"] # population size
			Q = pard["Q"] # institution size
			q = pard["q"] # strictness

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

			u_s = pard["u_s"] # mutation rate
			u_p = pard["u_p"] # performance error (e1)
			u_a = pard["u_a"] # assessment error (e2)

			num_gens = pard["num_gens"] # number of generations

			# name of output file
            output = pard["output"]

            println("--- running ", pard["nrun"], " --- ")
            flush(stdout)

			# initialize game and population
			game = Game(b, c, w, u_s, u_p, u_a, "pc")
			pop = institution_population(N, Q, q, game, reputation_norm, permitted_strategies)

			total_interactions = N^2

			coop_freq = Float64[]
			reputations = Float64[]
			strat_freqs = Array{Float64, 1}[]

			for g in 1:num_gens
				evolve!(pop)
				# allow num_gens/2 burn-in time, then start recording
				if g > num_gens/2
					# record coop frequency, reputations, and strategy frequencies
					push!(coop_freq, sum(pop.prev_actions)/total_interactions)
					push!(reputations, get_pub_reputations(pop))
					push!(strat_freqs, get_freqs(pop))
				end
			end

			# take averages
			coop_freq = mean(coop_freq)
			reputations = mean(reputations)
		    strat_freqs = [mean(hcat(strat_freqs...)[x,:]) for x in 1:4]

			pard["coop_freqs"] = coop_freq
			pard["reputations"] = reputations
			pard["strat_freqs"] = strat_freqs

            # return data to master process
            put!(results, pard)
        end
    end

    total_time_start = now()

    # load parameter sets into inputs channel
    nruns = 0
    for parset in parsets
        pard = Dict(zip(keys(pars), parset))
        #println(pard)
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
                 ["rep", "coop_freqs", "reputations", "strat_freqs", "seed1", "seed2", "seed3", "seed4"]...)
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
main(["--input", "submit/paper_institutions_noRDisc.json"])
