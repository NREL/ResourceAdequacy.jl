struct NonSequentialTemporalResultAccumulator{V,S,ES,SS} <: ResultAccumulator{V,S,ES,SS}
    droppedcount::Vector{SumVariance{V}}
    droppedsum::Vector{SumVariance{V}}
    system::S
    extractionspec::ES
    simulationspec::SS
    rngs::Vector{MersenneTwister}
end

function accumulator(extractionspec::ExtractionSpec,
                     simulationspec::SimulationSpec{NonSequential},
                     resultspec::Temporal, sys::SystemModel{N,L,T,P,E,V},
                     seed::UInt) where {N,L,T,P,E,V}

    nthreads = Threads.nthreads()
    nperiods = length(sys.timestamps)

    droppedcount = Vector{SumVariance{V}}(nperiods)
    droppedsum = Vector{SumVariance{V}}(nperiods)

    for t in 1:nperiods
        droppedcount[t] = Series(Sum(), Variance())
        droppedsum[t] = Series(Sum(), Variance())
    end

    rngs = Vector{MersenneTwister}(nthreads)
    rngs_temp = randjump(MersenneTwister(seed), nthreads)

    Threads.@threads for i in 1:nthreads
        rngs[i] = copy(rngs_temp[i])
    end

    return NonSequentialTemporalResultAccumulator(
        droppedcount, droppedsum,
        sys, extractionspec, simulationspec, rngs)

end

function update!(acc::NonSequentialTemporalResultAccumulator,
                 result::SystemOutputStateSummary, t::Int)

    fit!(acc.droppedcount[t], result.lolp_system)
    fit!(acc.droppedsum[t], sum(result.eue_regions))
    return

end

function update!(acc::NonSequentialTemporalResultAccumulator{V,SystemModel{N,L,T,P,E,V}},
                 sample::SystemOutputStateSample, t::Int, i::Int) where {N,L,T,P,E,V}

    shortfall = droppedload(sample)
    isshortfall = !isapprox(shortfall, 0.)
    droppedenergy = powertoenergy(shortfall, L, T, P, E)

    fit!(acc.droppedcount[t], V(isshortfall))
    fit!(acc.droppedsum[t], droppedenergy)

    return

end

function finalize(acc::NonSequentialTemporalResultAccumulator{V,<:SystemModel{N,L,T,P,E,V}}
                  ) where {N,L,T,P,E,V}

    timestamps = acc.system.timestamps
    nperiods = length(timestamps)

    if ismontecarlo(acc.simulationspec)

        # Accumulator summed results nsamples times, need to scale back down
        nsamples = acc.simulationspec.nsamples
        lolps =
            map(r -> LOLP{L,T}(r...), mean_stderr.(acc.droppedcount, nsamples))
        eues =
            map(r -> EUE{1,L,T,E}(r...), mean_stderr.(acc.droppedsum, nsamples))

    else

        # Accumulator summed once per timestep, no scaling required
        lolps = map(r -> LOLP{L,T}(r...), mean_stderr.(acc.droppedcount))
        eues = map(r -> EUE{1,L,T,E}(r...), mean_stderr.(acc.droppedsum))

    end

    return TemporalResult(timestamps, LOLE(lolps), lolps, EUE(eues), eues,
                          acc.extractionspec, acc.simulationspec)

end
