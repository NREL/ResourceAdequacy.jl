struct Spatial <: ResultSpec end

struct SpatialResultAccumulator{V,S,ES,SS} <: ResultAccumulator{V,S,ES,SS}
    droppedcount_overall::Vector{SumVariance{V}}
    droppedsum_overall::Vector{SumVariance{V}}
    droppedcount_region::Matrix{SumVariance{V}}
    droppedsum_region::Matrix{SumVariance{V}}
    localidx::Vector{Int}
    droppedcount_local::Vector{V}
    droppedsum_local::Vector{V}
    system::S
    extractionspec::ES
    simulationspec::SS
    rngs::Vector{MersenneTwister}
end

struct SpatialResult{
    N, # Number of timesteps simulated
    L, # Length of each timestep
    T <: Period, # Units of timestep duration
    E <: EnergyUnit, # Units for energy results
    V <: Real, # Numerical type of value data
    ES <: ExtractionSpec,
    SS <: SimulationSpec
} <: Result{N,L,T,V,ES,SS}

    regions::Vector{String}
    lole::LOLE{N,L,T,V}
    loles::Vector{LOLE{N,L,T,V}}
    eue::EUE{N,L,T,E,V}
    eues::Vector{EUE{N,L,T,E,V}}
    extractionspec::ES
    simulationspec::SS

    SpatialResult{}(
        regions::Vector{String},
        lole::LOLE{N,L,T,V}, loles::Vector{LOLE{N,L,T,V}},
        eue::EUE{N,L,T,E,V}, eues::Vector{EUE{N,L,T,E,V}},
        extractionspec::ES, simulationspec::SS) where {N,L,T,E,V,ES,SS} =
        new{N,L,T,E,V,ES,SS}(regions, lole, loles, eue, eues,
                             extractionspec, simulationspec)

end

LOLE(x::SpatialResult) = x.lole
LOLE(x::SpatialResult, r::Int) = x.loles[r]
LOLE(x::SpatialResult, r::AbstractString) =
    x.loles[findfirstunique(x.regions, r)]

EUE(x::SpatialResult) = x.eue
EUE(x::SpatialResult, r::Int) = x.eues[r]
EUE(x::SpatialResult, r::AbstractString) =
    x.eues[findfirstunique(x.regions, r)]

function accumulator(extractionspec::ExtractionSpec,
                     simulationspec::SimulationSpec,
                     resultspec::Spatial, sys::SystemModel{N,L,T,P,E,V},
                     seed::UInt) where {N,L,T,P,E,V}

    nthreads = Threads.nthreads()
    nregions = length(sys.regions)

    droppedcount_overall = Vector{SumVariance{V}}(nthreads)
    droppedsum_overall = Vector{SumVariance{V}}(nthreads)
    droppedcount_region = Matrix{SumVariance{V}}(nregions, nthreads)
    droppedsum_region = Matrix{SumVariance{V}}(nregions, nthreads)

    rngs = Vector{MersenneTwister}(nthreads)
    rngs_temp = randjump(MersenneTwister(seed), nthreads)
    localidx = zeros(Int, nthreads)
    localcount = Vector{V}(nthreads)
    localsum = Vector{V}(nthreads)

    Threads.@threads for i in 1:nthreads
        droppedcount_overall[i] = Series(Sum(), Variance())
        droppedsum_overall[i] = Series(Sum(), Variance())
        for r in 1:nregions
            droppedsum_region[r, i] = Series(Sum(), Variance())
            droppedcount_region[r, i] = Series(Sum(), Variance())
        end
        rngs[i] = copy(rngs_temp[i])
    end

    return SpatialResultAccumulator(
        droppedcount_overall, droppedsum_overall,
        droppedcount_region, droppedsum_region,
        localidx, localcount, localsum,
        sys, extractionspec, simulationspec, rngs)

end

function update!(acc::SpatialResultAccumulator{V},
                 sample::SystemOutputStateSample, t::Int, i::Int) where {V}

    error("Not yet implemented")
    return

end

function update!(acc::SpatialResultAccumulator,
                 result::SystemOutputStateSummary, t::Int)

    issequential(acc.simulationspec) &&
        error("Sequential analytical solutions are not currently supported.")

    thread = Threads.threadid()

    fit!(acc.droppedcount_overall[thread], result.lolp_system)
    fit!(acc.droppedsum_overall[thread], sum(result.eue_regions))

    for r in 1:length(acc.system.regions)
        fit!(acc.droppedcount_region[r, thread], result.lolp_regions[r])
        fit!(acc.droppedsum_region[r, thread], result.eue_regions[r])
    end

    return

end

function finalize(acc::SpatialResultAccumulator{V,<:SystemModel{N,L,T,P,E,V}}
                  ) where {N,L,T,P,E,V}

    regions = acc.system.regions
    nregions = length(regions)

    # Merge thread-local stats into final stats
    for i in 2:Threads.nthreads()

        merge!(acc.droppedcount_overall[1], acc.droppedcount[i])
        merge!(acc.droppedsum_overall[1], acc.droppedsum[i])

        for r in 1:nregions
            merge!(acc.droppedcount_region[r, 1], acc.droppedcount[r, i])
            merge!(acc.droppedsum_region[r, 1], acc.droppedsum[r, i])
        end

    end

    if ismontecarlo(acc.simulationspec)

        # Accumulator summed results nsamples times, need to scale back down
        nsamples = acc.simulationspec.nsamples
        lole = LOLE{N,L,T}(mean_stderr(acc.droppedcount_overall[1], nsamples)...)
        loles = map(r -> LOLE{N,L,T}(r...),
                    mean_stderr.(acc.droppedcount_region[:, 1], nsamples))
        eue = EUE{N,L,T,E}(mean_stderr(acc.droppedsum_overall[1], nsamples)...)
        eues = map(r -> EUE{N,L,T,E}(r...),
                   mean_stderr.(acc.droppedsum_region[:, 1], nsamples))

    else

        # Accumulator summed once per timestep, no scaling required
        lole = LOLE{N,L,T}(mean_stderr(acc.droppedcount_overall[1])...)
        loles = map(r -> LOLE{N,L,T}(r...),
                    mean_stderr.(acc.droppedcount_region[:, 1]))
        eue = EUE{N,L,T,E}(mean_stderr(acc.droppedsum_overall[1])...)
        eues = map(r -> EUE{N,L,T,E}(r...),
                   mean_stderr.(acc.droppedsum_region[:, 1]))

    end

    return SpatialResult(regions, lole, loles, eue, eues,
                          acc.extractionspec, acc.simulationspec)

end
