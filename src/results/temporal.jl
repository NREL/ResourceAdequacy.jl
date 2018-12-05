struct Temporal <: ResultSpec end

struct TemporalResultAccumulator{V,S,ES,SS} <: ResultAccumulator{V,S,ES,SS}
    droppedcount_overall::Vector{SumVariance{V}}
    droppedsum_overall::Vector{SumVariance{V}}
    droppedcount_period::Matrix{SumVariance{V}}
    droppedsum_period::Matrix{SumVariance{V}}
    localidx::Vector{Int}
    droppedcount_local::Vector{V}
    droppedsum_local::Vector{V}
    system::S
    extractionspec::ES
    simulationspec::SS
    rngs::Vector{MersenneTwister}
end

struct TemporalResult{
    N, # Number of timesteps simulated
    L, # Length of each timestep
    T <: Period, # Units of timestep duration
    E <: EnergyUnit, # Units for energy results
    V <: Real, # Numerical type of value data
    ES <: ExtractionSpec,
    SS <: SimulationSpec
} <: Result{N,L,T,V,ES,SS}

    timestamps::StepRange{DateTime,T}
    lole::LOLE{N,L,T,V}
    lolps::Vector{LOLP{L,T,V}}
    eue::EUE{N,L,T,E,V}
    eues::Vector{EUE{1,L,T,E,V}}
    extractionspec::ES
    simulationspec::SS

    TemporalResult{}(
        timestamps::StepRange{DateTime,T},
        lole::LOLE{N,L,T,V}, lolps::Vector{LOLP{L,T,V}},
        eue::EUE{N,L,T,E,V}, eues::Vector{EUE{1,L,T,E,V}},
        extractionspec::ES, simulationspec::SS) where {N,L,T,E,V,ES,SS} =
        new{N,L,T,E,V,ES,SS}(timestamps, lole, lolps, eue, eues,
                             extractionspec, simulationspec)

end

LOLE(x::TemporalResult) = x.lole
LOLP(x::TemporalResult, t::Int) = x.lolps[t]
LOLP(x::TemporalResult, dt::DateTime) =
    x.lolps[findfirstunique(x.timestamps, dt)]

EUE(x::TemporalResult) = x.eue
EUE(x::TemporalResult, t::Int) = x.eues[t]
EUE(x::TemporalResult, dt::DateTime) =
    x.eues[findfirstunique(x.timestamps, dt)]

function accumulator(extractionspec::ExtractionSpec,
                     simulationspec::SimulationSpec,
                     resultspec::Temporal, sys::SystemModel{N,L,T,P,E,V},
                     seed::UInt) where {N,L,T,P,E,V}

    nthreads = Threads.nthreads()
    nperiods = length(sys.timestamps)

    droppedcount_overall = Vector{SumVariance{V}}(nthreads)
    droppedsum_overall = Vector{SumVariance{V}}(nthreads)
    droppedcount_period = Matrix{SumVariance{V}}(nperiods, nthreads)
    droppedsum_period = Matrix{SumVariance{V}}(nperiods, nthreads)

    rngs = Vector{MersenneTwister}(nthreads)
    rngs_temp = randjump(MersenneTwister(seed), nthreads)
    localidx = zeros(Int, nthreads)
    localcount = Vector{V}(nthreads)
    localsum = Vector{V}(nthreads)

    Threads.@threads for i in 1:nthreads
        droppedcount_overall[i] = Series(Sum(), Variance())
        droppedsum_overall[i] = Series(Sum(), Variance())
        for t in 1:nperiods
            droppedsum_period[t, i] = Series(Sum(), Variance())
            droppedcount_period[t, i] = Series(Sum(), Variance())
        end
        rngs[i] = copy(rngs_temp[i])
    end

    return TemporalResultAccumulator(
        droppedcount_overall, droppedsum_overall,
        droppedcount_period, droppedsum_period,
        localidx, localcount, localsum,
        sys, extractionspec, simulationspec, rngs)

end

function update!(acc::TemporalResultAccumulator{V},
                 sample::SystemOutputStateSample, t::Int, i::Int) where {V}

    error("Not yet implemented")
    return

end

function update!(acc::TemporalResultAccumulator,
                 result::SystemOutputStateSummary, t::Int)

    issequential(acc.simulationspec) &&
        error("Sequential analytical solutions are not currently supported.")

    thread = Threads.threadid()

    fit!(acc.droppedcount_overall[thread], result.lolp_system)
    fit!(acc.droppedsum_overall[thread], sum(result.eue_regions))

    fit!(acc.droppedcount_period[t, thread], result.lolp_system)
    fit!(acc.droppedsum_period[t, thread], sum(result.eue_regions))

    return

end

function finalize(acc::TemporalResultAccumulator{V,<:SystemModel{N,L,T,P,E,V}}
                  ) where {N,L,T,P,E,V}

    timestamps = acc.system.timestamps
    nperiods = length(timestamps)

    # Merge thread-local stats into final stats
    for i in 2:Threads.nthreads()

        merge!(acc.droppedcount_overall[1], acc.droppedcount[i])
        merge!(acc.droppedsum_overall[1], acc.droppedsum[i])

        for t in 1:nperiods
            merge!(acc.droppedcount_period[t, 1], acc.droppedcount[t,i])
            merge!(acc.droppedsum_period[t, 1], acc.droppedsum[t,i])
        end

    end

    if ismontecarlo(acc.simulationspec)

        # Accumulator summed results nsamples times, need to scale back down
        nsamples = acc.simulationspec.nsamples
        lole = LOLE{N,L,T}(mean_stderr(acc.droppedcount_overall[1], nsamples)...)
        lolps = map((m,s) -> LOLP{L,T}(m,s),
                    mean_stderr.(acc.droppedcount_period[:, 1], nsamples))
        eue = EUE{N,L,T,E}(mean_stderr(acc.droppedsum_overall[1], nsamples)...)
        eues = map((m,s) -> EUE{1,L,T,E}(m,s),
                   mean_stderr.(acc.droppedsum_period[:, 1], nsamples))

    else

        # Accumulator summed once per timestep, no scaling required
        lole = LOLE{N,L,T}(mean_stderr(acc.droppedcount_overall[1])...)
        lolps = map(r -> LOLP{L,T}(r[1],r[2]),
                    mean_stderr.(acc.droppedcount_period[:, 1]))
        eue = EUE{N,L,T,E}(mean_stderr(acc.droppedsum_overall[1])...)
        eues = map(r -> EUE{1,L,T,E}(r[1], r[2]),
                   mean_stderr.(acc.droppedsum_period[:, 1]))

    end

    return TemporalResult(timestamps, lole, lolps, eue, eues,
                          acc.extractionspec, acc.simulationspec)

end
