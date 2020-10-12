#!/usr/bin/env julia

## parallel_lending.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Parallelized implementation of Institutions.

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

    s = ArgParseSettings(description =
        "run ReputationSets simulations across multiple cores")
    @add_arg_table s begin
        "--ncpus"
            arg_type = Int64
            default = max(round(Int, Sys.CPU_THREADS / 2), 1)
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
		"u_s"     => Dict("value" => 0.01,     "type" => Float64),
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

            N = pard["N"]
			E = pard["E"]
			Q = pard["Q"]
			q = pard["q"]

			reputation_norm = pard["reputation_norm"]

			if pard["permitted_strategies"] == "all"
				permitted_strategies = [1,2,3,4]
			else
				permitted_strategies = parse.(Int64, split(pard["permitted_strategies"], ""))
			end

			b = pard["b"]
			c = pard["c"]
			w = pard["w"]
			independent_board = pard["independent_board"]


			u_s = pard["u_s"]
			u_p = pard["u_p"]
			u_a = pard["u_a"]

			run_type = pard["run_type"]

			num_trials = pard["num_trials"]
			runs_per_trial = pard["runs_per_trial"]

            output = pard["output"]

            println("--- running ", pard["nrun"], " --- ")
            flush(stdout)

			game = Game(b, c, w, u_s, u_p, u_a, "pc")

			total_interactions = N^2

			successes = 0
			failures = 0

			for g in 1:runs_per_trial
				pop = fixation_population(N, E, Q, q, game, reputation_norm, permitted_strategies, independent_board)
				if run_type == "equal"
					pop.strategies[1:16] .= 1
					pop.strategies[17:33] .= 2
					pop.strategies[34:50] .= 3
				elseif run_type == "AllC"
					pop.strategies[1:46] .= 1
					pop.strategies[47:48] .= 2
					pop.strategies[49:50] .= 3
				elseif run_type == "AllD"
					pop.strategies[1:2] .= 1
					pop.strategies[3:48] .= 2
					pop.strategies[49:50] .= 3
				elseif run_type == "Disc"
					pop.strategies[1:2] .= 1
					pop.strategies[3:4] .= 2
					pop.strategies[5:50] .= 3
				else
				end
				if sum([pop.strategies[x] .== 3 for x in 1:pop.N]) > 0
					pop.is_follower[rand(filter(x -> pop.strategies[x] .== 3, 1:pop.N))] = 1
				else
					x = rand(1:pop.N)
					pop.strategies[x] = 3
					pop.is_follower[x] = 1
				end
				while sum(pop.is_follower) ∉ [0, pop.N]
					evolve!(pop)
				end
				if sum(pop.is_follower) == pop.N
					successes += 1
				else
					failures += 1
				end
			end


			pard["successes"] = successes
			pard["failures"] = failures

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
                 ["rep", "successes", "failures", "seed1", "seed2", "seed3", "seed4"]...)
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

#main(ARGS)

main(["--input", "submit/fixation_paper_var_b.json"])