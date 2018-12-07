struct NonSequentialSpatialResultAccumulator{V,S,ES,SS} <: ResultAccumulator{V,S,ES,SS}

    droppedcount_overall_valsum::Vector{V}
    droppedcount_overall_varsum::Vector{V}
    droppedsum_overall_valsum::Vector{V}
    droppedsum_overall_varsum::Vector{V}
    droppedcount_region_valsum::Matrix{V}
    droppedcount_region_varsum::Matrix{V}
    droppedsum_region_valsum::Matrix{V}
    droppedsum_region_varsum::Matrix{V}

    periodidx::Vector{Int}
    droppedcount_overall_period::Vector{SumVariance{V}}
    droppedsum_overall_period::Vector{SumVariance{V}}
    droppedcount_region_period::Matrix{SumVariance{V}}
    droppedsum_region_period::Matrix{SumVariance{V}}

    system::S
    extractionspec::ES
    simulationspec::SS
    rngs::Vector{MersenneTwister}

end

function accumulator(extractionspec::ExtractionSpec,
                     simulationspec::SimulationSpec{NonSequential},
                     resultspec::Spatial, sys::SystemModel{N,L,T,P,E,V},
                     seed::UInt) where {N,L,T,P,E,V}

    nthreads = Threads.nthreads()
    nregions = length(sys.regions)

    droppedcount_overall_valsum = zeros(nthreads)
    droppedcount_overall_varsum = zeros(nthreads)
    droppedsum_overall_valsum = zeros(nthreads)
    droppedsum_overall_varsum = zeros(nthreads)

    droppedcount_region_valsum = zeros(nregions, nthreads)
    droppedcount_region_varsum = zeros(nregions, nthreads)
    droppedsum_region_valsum = zeros(nregions, nthreads)
    droppedsum_region_varsum = zeros(nregions, nthreads)

    periodidx = zeros(Int, nthreads)
    droppedcount_overall_period = Vector{SumVariance{V}}(nthreads)
    droppedsum_overall_period = Vector{SumVariance{V}}(nthreads)
    droppedcount_region_period = Matrix{SumVariance{V}}(nregions, nthreads)
    droppedsum_region_period = Matrix{SumVariance{V}}(nregions, nthreads)

    rngs = Vector{MersenneTwister}(nthreads)
    rngs_temp = randjump(MersenneTwister(seed), nthreads)

    Threads.@threads for i in 1:nthreads
        droppedcount_overall_period[i] = Series(Sum(), Variance())
        droppedsum_overall_period[i] = Series(Sum(), Variance())
        for r in 1:nregions
            droppedcount_region_period[r, i] = Series(Sum(), Variance())
            droppedsum_region_period[r, i] = Series(Sum(), Variance())
        end
        rngs[i] = copy(rngs_temp[i])
    end

    return NonSequentialSpatialResultAccumulator(
        droppedcount_overall_valsum, droppedcount_overall_varsum,
        droppedsum_overall_valsum, droppedsum_overall_varsum,
        droppedcount_region_valsum, droppedcount_region_varsum,
        droppedsum_region_valsum, droppedsum_region_varsum,
        periodidx, droppedcount_overall_period, droppedsum_overall_period,
        droppedcount_region_period, droppedsum_region_period,
        sys, extractionspec, simulationspec, rngs)

end

function update!(acc::NonSequentialSpatialResultAccumulator,
                 result::SystemOutputStateSummary, t::Int)

    thread = Threads.threadid()

    acc.droppedcount_overall_valsum[thread] += result.lolp_system
    acc.droppedsum_overall_valsum[thread] += sum(result.eue_regions)

    for r in 1:length(acc.system.regions)
        acc.droppedcount_region_valsum[r, thread] += result.lolp_regions[r]
        acc.droppedsum_region_valsum[r, thread] += result.eue_regions[r]
    end

    return

end

function update!(acc::NonSequentialSpatialResultAccumulator{V,SystemModel{N,L,T,P,E,V}},
                 sample::SystemOutputStateSample, t::Int, i::Int) where {N,L,T,P,E,V}

    thread = Threads.threadid()
    nregions = length(acc.system.regions)

    if t != acc.periodidx[thread]

        # Previous local period has finished,
        # so store previous local result and reset

        transferperiodresults!(
            acc.droppedcount_overall_valsum, acc.droppedcount_overall_varsum,
            acc.droppedcount_overall_period, thread)
        
        # periodcount_sum, periodcount_var =
            # value(acc.droppedcount_overall_local[thread])
        # acc.droppedcount_overall_valsum[thread] += periodcount_sum
        # acc.droppedcount_overall_varsum[thread] += periodcount_var
        # acc.droppedcount_overall_local[thread] = Series(Sum(), Variance())

        transferperiodresults!(
            acc.droppedsum_overall_valsum, acc.droppedsum_overall_varsum,
            acc.droppedsum_overall_period, thread)

        # periodsum_sum, periodsum_var =
        #     value(acc.droppedsum_overall_local[thread])
        # acc.droppedsum_overall_valsum[thread] += periodsum_sum
        # acc.droppedsum_overall_varsum[thread] += periodsum_var
        # acc.droppedsum__overall_local[thread] = Series(Sum(), Variance())

        for r in 1:length(acc.system.regions)

            transferperiodresults!(
                acc.droppedcount_region_valsum, acc.droppedcount_region_varsum,
                acc.droppedcount_region_period, r, thread)

            # periodcount_region_sum, periodcount_region_var =
            #     value(acc.droppedcount_region_local[r, thread])
            # acc.droppedcount_region_valsum[r, thread] += periodcount_region_sum
            # acc.droppedcount_region_varsum[r, thread] += periodcount_region_var
            # acc.droppedcount_region_local[r, thread] = Series(Sum(), Variance())

            transferperiodresults!(
                acc.droppedsum_region_valsum, acc.droppedsum_region_varsum,
                acc.droppedsum_region_period, r, thread)

            # periodsum_region_sum, periodsum_region_var =
            #     value(acc.droppedsum_region_local[r, thread])
            # acc.droppedsum_region_valsum[r, thread] += periodsum_region_sum
            # acc.droppedsum_region_varsum[r, thread] += periodsum_region_var
            # acc.droppedsum_region_local[r, thread] = Series(Sum(), Variance())

        end

        acc.periodidx[thread] = t

    end

    shortfalls = droppedloads(sample)
    shortfall = sum(shortfalls)

    fit!(acc.droppedcount_overall_period[thread], approxnonzero(shortfall))
    fit!(acc.droppedsum_overall_period[thread],
         powertoenergy(shortfall, L, T, P, E))

    for r in 1:nregions
        fit!(acc.droppedcount_region_period[r, thread],
             approxnonzero(shortfalls[r]))
        fit!(acc.droppedsum_region_period[r, thread],
             powertoenergy(shortfalls[r], L, T, P, E))
    end

    return

end

function finalize(acc::NonSequentialSpatialResultAccumulator{V,<:SystemModel{N,L,T,P,E,V}}
                  ) where {N,L,T,P,E,V}

    regions = acc.system.regions
    nregions = length(regions)

    # Add the final local results
    for thread in 1:Threads.threadid()

        transferperiodresults!(
            acc.droppedcount_overall_valsum, acc.droppedcount_overall_varsum,
            acc.droppedcount_overall_period, thread)

        transferperiodresults!(
            acc.droppedsum_overall_valsum, acc.droppedsum_overall_varsum,
            acc.droppedsum_overall_period, thread)

    end

    if ismontecarlo(acc.simulationspec)

        # Accumulator summed results nsamples times, need to scale back down
        nsamples = acc.simulationspec.nsamples

        lole = sum(acc.droppedcount_overall_valsum) / nsamples
        lole_stderr = sqrt(sum(acc.droppedcount_overall_varsum)) / nsamples

        loles = vec(sum(acc.droppedcount_region_valsum, 2)) ./ nsamples
        loles_stderr = sqrt.(vec(sum(acc.droppedcount_region_varsum, 2))) ./ nsamples

        eue = sum(acc.droppedsum_overall_valsum) / nsamples
        eue_stderr = sqrt(sum(acc.droppedsum_overall_varsum)) / nsamples 

        eues = vec(sum(acc.droppedsum_region_valsum, 2)) ./ nsamples
        eues_stderr = sqrt.(vec(sum(acc.droppedsum_region_varsum, 2))) ./ nsamples 

    else

        # Accumulator summed once per timestep, no scaling required

        lole = sum(acc.droppedcount_overall_valsum)
        loles = vec(sum(acc.droppedcount_overall_valsum, 2))
        lole_stderr = loles_stderr = zero(V)

        eue = sum(acc.droppedsum_overall_valsum)
        eues = vec(sum(acc.droppedsum_overall_valsum, 2))
        eue_stderr = eues_stderr = zero(V)

    end

    return SpatialResult(regions,
                         LOLE{N,L,T}(lole, lole_stderr),
                         LOLE{N,L,T}.(loles, loles_stderr),
                         EUE{N,L,T,E}(eue, eue_stderr),
                         EUE{N,L,T,E}.(eues, eues_stderr),
                         acc.extractionspec, acc.simulationspec)

end
