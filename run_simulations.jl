@everywhere using MathProgBase
@everywhere using Clp
@everywhere using RowEchelon
@everywhere using DataFrames
@everywhere using CSV
@everywhere using JLD
@everywhere include("sequentialsinglenode.jl")

paramslist = [SimulationParams(0.1, 0.1, drlevel/100, 1000, "results_DR$(drlevel).jld")
              for drlevel in 0:25:100]
pmap(runsimulation, paramslist)
