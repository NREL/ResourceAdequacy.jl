struct Minimal <: ResultSpec end

# TODO: Need to enforce consistency between V and SystemModel{.., V}
struct MinimalResultAccumulator{V,S,SS} <: ResultAccumulator{V,S,SS}
    droppedcount::Vector{SumVariance{V}}
    droppedsum::Vector{SumVariance{V}}
    system::S
    simulationspec::SS
    rngs::Vector{MersenneTwister}
end

struct MinimalResult{
    N, # Number of timesteps simulated
    L, # Length of each timestep
    T <: Period, # Units of timestep duration
    E <: EnergyUnit, # Units for energy results
    V <: Real, # Numerical type of value data
    ES <: ExtractionSpec,
    SS <: SimulationSpec
} <: Result{N,L,T,V,ES,SS}

    lole::LOLE{N,L,T,V}
    eue::EUE{E,N,L,T,V}
    extractionspec::ES
    simulationspec::SS

    MinimalResult{}(
        lole::LOLE{N,L,T,V}, eue::EUE{E,N,L,T,V},
        extractionspec::ES, simulationspec::SS) where {N,L,T,E,V,ES,SS} =
        new{N,L,T,E,V,ES,SS}(lole, eue, extractionspec, simulationspec)

end

LOLE(x::MinimalResult) = x.lole
EUE(x::MinimalResult) = x.eue

function accumulator(extractionspec::ExtractionSpec,
                     simulationspec::SimulationSpec,
                     resultspec::Minimal, sys::SystemModel{N,L,T,P,E,V},
                     seed::UInt) where {N,L,T,P,E,V}

    nthreads = Threads.nthreads()

    droppedcount = Vector{SumVariance{V}}(n_threads)
    droppedsum = Vector{SumVariance{V}}(n_threads)
    rngs = Vector{MersenneTwister}(n_threads)
    rngs_temp = randjump(MersenneTwister(seed), n_threads)

    Threads.@threads for i in 1:nthreads
        droppedcount[i] = Series(Sum(), Variance())
        droppedsum[i] = Series(Sum(), Variance())
        rngs[i] = copy(rngs_temp[i])
    end

    return MinimalResultAccumulator(droppedcount, droppedsum, sys,
                                    simulationspec, rngs)

end

function update!(acc::MinimalResultAccumulator{V},
                 unservedenergy::V,
                 unservedenergyperiods::V) where {V<:Real}

    issequential(acc.simulationspec) || return

    fit!(acc.droppedsum, unservedenergy)
    fit!(acc.droppedcount, unservedenergyperiods)
    return

end

function update!(acc::MinimalResultAccumulator{V},
                 unservedenergy::V,
                 unservedenergyperiods::V, ::Int) where {V<:Real}

    issequential(acc.simulationspec) && return

    fit!(acc.droppedsum, unservedenergy)
    fit!(acc.droppedcount, unservedenergyperiods)
    return

end

function finalize(acc::MinimalResultAccumulator{V,<:SystemModel{N,L,T,P,E,V}},
                  extractionspec::ExtractionSpec) where {N,L,T,P,E,V}

    # TODO: Merge thread stats into final stats

    if ismontecarlo(acc.simulationspec)

        nsamples = acc.simulationspec.nsamples
        lole, lole_stderr = mean_stderr(acc.droppedcount, nsamples)
        eue, eue_stderr = mean_stderr(acc.droppedsum, nsamples)

    else

        lole, lole_stderr = mean_stderr(acc.droppedcount)
        eue, eue_stderr = mean_stderr(acc.droppedsum)

    end
    
    return MinimalResult{P}(
        LOLE{N,L,T}(lole, lole_stderr),
        EUE{E,N,L,T}(eue, eue_stderr),
        extractionspec, acc.simulationspec)

end
