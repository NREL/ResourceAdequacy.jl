@everywhere using MathProgBase
@everywhere using Clp
@everywhere using RowEchelon
@everywhere using DataFrames
@everywhere using CSV
@everywhere using JLD

@everywhere include("SimulationParams.jl")
@everywhere include("ResultStore.jl")
@everywhere include("sequentialsinglenode.jl")
@everywhere include("capacityvalue.jl")

mode = ARGS[1]

if mode == "RA"

	# Run resource adequacy assessments
	paramslist = [SimulationParams(0.1, 0.1, drlevel/100, 0.2, 1000, "results_DR$(drlevel).jld", 0.)
		      for drlevel in 0:25:100]
	pmap(runsimulation, paramslist)

elseif mode == "CV"

	# Run capacity valuations
	targetloles = [0.06, 0.14, 0.33, 1.11]
	capacityvalues = capacityvalue(:lole, targetloles, 0.01, 3602., 22, "cvcurve_lole.csv")

	println("\nEquivalent Firm Capacities:")
	foreach((drfrac, cv) -> println(drfrac/100, " => ", cv, " MW"),
		25:25:100,  capacityvalues)

else

    warn("Unrecognized command, should be either " *
         "RA (resource adequacy assessment) or CV (capacity valuation")

end
