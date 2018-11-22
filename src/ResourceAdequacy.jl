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

# Abstract spec interfaces and instances
spec_instances = [
    ("extraction", ["backcast", "repra"]),
    ("simulation", ["nonsequentialcopperplate", "nonsequentialnetworkflow"]),
    ("result", ["minimal", "temporal", "spatial", "spatiotemporal", "network"])
]

# Load abstract interfaces
for (spec, _) in specs
    include("abstractspecs/" * spec * ".jl")
end

# Load concrete instances
for (spec, instances) in specs, instance in instances
    include(spec * "s/" * instance * ".jl")
end

end # module
