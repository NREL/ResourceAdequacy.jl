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

# Basic functionality
include("utils/utils.jl")
include("systemdata/systemdata.jl")
include("metrics/metrics.jl")
include("abstractspecs/abstractspecs.jl")

# Spec instances
spec_instances = [
    ("extraction", ["backcast", "repra"]),
    ("simulation", ["nonsequentialcopperplate"]), #, "nonsequentialnetworkflow"]),
    ("result", ["minimal"])  #, "temporal", "spatial", "spatiotemporal", "network"])
]
for (spec, instances) in spec_instances, instance in instances
    include(spec * "s/" * instance * ".jl")
end

end # module
