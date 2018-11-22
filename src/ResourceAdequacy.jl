module ResourceAdequacy

using Base.Dates
using StatsBase
using OnlineStats
using Distributions
using LightGraphs
using MathProgBase

export

    assess,

    # Units
    Year, Month, Week, Day, Hour, Minute,
    MW, GW,
    MWh, GWh, TWh,

    # Metrics
    LOLP, LOLE, EUE,
    val, stderr,

    # Distribution extraction specifications
    Backcast, REPRA,

    # Simulation specifications
    NonSequentialCopperplate, SequentialCopperplate,
    NonSequentialNetworkFlow, SequentialNetworkFlow,

    # Result specifications
    MinimalResult, NetworkResult,

    # Result methods
    timestamps

CapacityDistribution{T} = Distributions.Generic{T,Float64,Vector{T}}
CapacitySampler{T} = Distributions.GenericSampler{T, Vector{T}}

include("utils.jl")
include("metrics.jl")
include("extraction.jl")
include("simulation.jl")
include("results.jl")
include("systemdata.jl")

end # module
