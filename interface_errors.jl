#libraries for loading the test system
using HDF5
using JLD
using Dates #this was key for avoiding the warning about datetime elements being missing (which caused things to fail subsequently)

#read in the test system
cd("/Users/mschwarz/") #filepath specific to my system
file = load("RTS_10.jld")
mysystemmodel = file["DAY_AHEAD"] #extract test system from file

##libraries for doing the RA analysis
using ResourceAdequacy
using PRAS

##The interface flow and utilization factor queries only works with certain combinations of modules. 
##As instructed in the "ResourceAdequacy.jl" README, these queries are only relevant when the "Network()"
##result spec is used.

##These analyses work:
works1 = assess(Backcast(), NonSequentialCopperplate(), Network(), mysystemmodel)
works2 = assess(REPRA(1,10), NonSequentialCopperplate(), Network(), mysystemmodel)
flow_1 = ExpectedInterfaceFlow(works1, "1", "2", DateTime(2020,4,27,13))
flow_2 = ExpectedInterfaceFlow(works2, "1", "2", DateTime(2020,4,27,13))
utilization1 = ExpectedInterfaceUtilization(works1, "1", "2", DateTime(2020,4,27,13))
utilization2 = ExpectedInterfaceUtilization(works2, "1", "2", DateTime(2020,4,27,13))

##These analyses don't work:
error1=assess(Backcast(), NonSequentialNetworkFlow(100), Network(), mysystemmodel)
error2=assess(REPRA(1, 10), NonSequentialNetworkFlow(100), Network(), mysystemmodel)

#So the issue seems to be with assess(), not ExpectedInterfaceFlow() or ExpectedInterfaceUtilization().

#They yield the same error message:
# ERROR: NaN is not a valid interface utilization expectation
# Stacktrace:
#  [1] error(::String) at ./error.jl:33
#  [2] ExpectedInterfaceUtilization{1,1,Hour,V} where V<:Real(::Float64, ::Float64) at /Users/mschwarz/.julia/packages/ResourceAdequacy/L0wVg/src/metrics/ExpectedInterfaceUtilization.jl:8
#  [3] makemetric(::Type, ::OnlineStatsBase.Series{Number,Tuple{OnlineStatsBase.Mean{Float64,OnlineStatsBase.EqualWeight},OnlineStatsBase.Variance{Float64,OnlineStatsBase.EqualWeight}}}) at /Users/mschwarz/.julia/packages/ResourceAdequacy/L0wVg/src/utils/misc.jl:20
#  [4] _broadcast_getindex_evalf at ./broadcast.jl:578 [inlined]
#  [5] _broadcast_getindex at ./broadcast.jl:561 [inlined]
#  [6] getindex at ./broadcast.jl:511 [inlined]
#  [7] macro expansion at ./broadcast.jl:843 [inlined]
#  [8] macro expansion at ./simdloop.jl:73 [inlined]
#  [9] copyto! at ./broadcast.jl:842 [inlined]
#  [10] copyto! at ./broadcast.jl:797 [inlined]
#  [11] copy at ./broadcast.jl:773 [inlined]
#  [12] materialize at ./broadcast.jl:753 [inlined]
#  [13] finalize(::ResourceAdequacy.NonSequentialNetworkResultAccumulator{Float64,PRASBase.SystemModel{8784,1,Hour,MW,MWh,Float64},REPRA,NonSequentialNetworkFlow}) at /Users/mschwarz/.julia/packages/ResourceAdequacy/L0wVg/src/results/network_nonsequentialaccumulator.jl:155
#  [14] assess(::REPRA, ::NonSequentialNetworkFlow, ::Network, ::PRASBase.SystemModel{8784,1,Hour,MW,MWh,Float64}, ::UInt64) at /Users/mschwarz/.julia/packages/ResourceAdequacy/L0wVg/src/abstractspecs/simulation.jl:54
#  [15] assess(::REPRA, ::NonSequentialNetworkFlow, ::Network, ::PRASBase.SystemModel{8784,1,Hour,MW,MWh,Float64}) at /Users/mschwarz/.julia/packages/ResourceAdequacy/L0wVg/src/abstractspecs/simulation.jl:43
#  [16] top-level scope at none:0


