struct NonSequentialNetworkFlow <: SimulationSpec{NonSequential}
    nsamples::Int

    function NonSequentialNetworkFlow(nsamples::Int)
        @assert nsamples > 0
        new(nsamples)
    end
end

iscopperplate(::NonSequentialNetworkFlow) = false

accumulator(ss::NonSequentialNetworkFlow, rs::MinimalResult,
            sys::SystemStateDistribution) =
                SinglePeriodNetworkMinimalResultAccumulator(ss, sys)

accumulator(ss::NonSequentialNetworkFlow, rs::NetworkResult,
            sys::SystemStateDistribution) =
                SinglePeriodNetworkResultAccumulator(ss, rs, sys)

function assess(simulationspec::NonSequentialNetworkFlow,
                resultspec::ResultSpec,
                system::SystemStateDistribution{N,T,P,E,Float64}) where {N,T,P,E}

    systemsampler = SystemSampler(system)
    sink_idx = nv(systemsampler.graph)
    source_idx = sink_idx-1
    n = sink_idx-2

    resultaccumulator = accumulator(simulationspec, resultspec, system)

    state_matrix = zeros(sink_idx, sink_idx)
    flow_matrix = Array{Float64}(sink_idx, sink_idx)
    height = Array{Int}(sink_idx)
    count = Array{Int}(2*sink_idx+1)
    excess = Array{Float64}(sink_idx)
    active = Array{Bool}(sink_idx)

    for i in 1:simulationspec.nsamples

        rand!(state_matrix, systemsampler)
        systemload, flow_matrix =
            LightGraphs.push_relabel!(
                flow_matrix, height, count, excess, active,
                systemsampler.graph, source_idx, sink_idx, state_matrix)

        savesample!(resultaccumulator, state_matrix, flow_matrix, sink_idx, n)

    end

    # return finalize(resultaccumulator)

end

function all_load_served(A::Matrix{T}, B::Matrix{T}, sink::Int, n::Int) where T
    served = true
    i = 1
    while served && (i <= n)
        served = A[i, sink] ≈ B[i, sink]
        i += 1
    end
    return served
end
