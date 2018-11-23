struct Minimal <: ResultSpec end

# TODO: Need to enforce consistency between V and SystemModel{.., V}
struct MinimalResultAccumulator{V,S,SS} <: ResultAccumulator{V,S,SS}
    droppedcount::Vector{SumVariance{V}}
    droppedsum::Vector{SumVariance{V}}
    localidx::Vector{Int}
    droppedcount_local::Vector{V}
    droppedsum_local::Vector{V}
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

    droppedcount = Vector{SumVariance{V}}(nthreads)
    droppedsum = Vector{SumVariance{V}}(nthreads)
    rngs = Vector{MersenneTwister}(nthreads)
    rngs_temp = randjump(MersenneTwister(seed), nthreads)

    localidx = zeros(Int, nthreads)
    localcount = Vector{Int}(nthreads)
    localsum = Vector{V}(nthreads)

    Threads.@threads for i in 1:nthreads
        droppedcount[i] = Series(Sum(), Variance())
        droppedsum[i] = Series(Sum(), Variance())
        rngs[i] = copy(rngs_temp[i])
    end

    return MinimalResultAccumulator(
        droppedcount, droppedsum, localidx, localcount, localsum,
        sys, simulationspec, rngs)

end

function update!(acc::MinimalResultAccumulator{V},
                 sample::SystemOutputStateSample, t::Int, i::Int) where {V}

    thread = Threads.threadid()
    droppedload = _(sample) #TODO: Find this function
    isshortfall = isapprox(droppedload, 0.)
    droppedenergy = powertoenergy(droppedload, N, T, P, E)

    # Sequential simulations are grouped by simulation,
    # while NonSequential are grouped by timestep
    localidx = issequential(acc.simulationspec) ? i : t

    if localidx != acc.localidx[thread]

        # Previous local simulation/timestep has finished,
        # so store previous local result and reset

        fit!(acc.droppedcount[thread], acc.droppedcount_local[thread])
        fit!(acc.droppedsum[thread], acc.droppedsum_local[thread])

        acc.localidx[thread] = i
        acc.droppedcount_local[thread] = V(isshortfall)
        acc.droppedsum_local[thread] = droppedenergy

    elseif isshortfall

        # Local simulation/timestep is still ongoing
        # Load was dropped, update local tracking

        acc.droppedcount_local[thread] += one(V)
        acc.droppedsum_local[thread] += droppedenergy

    end

    return

end

function update!(acc::MinimalResultAccumulator,
                 result::SystemOutputStateSummary, t::Int)

    issequential(acc.simulationspec) &&
        error("Sequential analytical solutions are not currently supported.")

    thread = Threads.threadid()
    fit!(acc.droppedsum[thread], result.lolp_system)
    fit!(acc.droppedcount[thread], sum(result.eue_regions))

    return

end

function finalize(acc::MinimalResultAccumulator{V,<:SystemModel{N,L,T,P,E,V}},
                  extractionspec::ExtractionSpec) where {N,L,T,P,E,V}

    # Merge thread-local stats into final stats
    for i in 2:Threads.nthreads()
        merge!(acc.droppedcount[1], acc.droppedcount[i])
        merge!(acc.droppedsum[1], acc.droppedsum[i])
    end

    if ismontecarlo(acc.simulationspec)
        # Accumulator summed results nsamples times, to scale back down
        nsamples = acc.simulationspec.nsamples
        lole, lole_stderr = mean_stderr(acc.droppedcount[1], nsamples)
        eue, eue_stderr = mean_stderr(acc.droppedsum[1], nsamples)
    else
        # Accumulator summed once per timestep, no scaling required
        lole, lole_stderr = mean_stderr(acc.droppedcount[1])
        eue, eue_stderr = mean_stderr(acc.droppedsum[1])
    end

    return MinimalResult(
        LOLE{N,L,T}(lole, lole_stderr),
        EUE{E,N,L,T}(eue, eue_stderr),
        extractionspec, acc.simulationspec)

end
