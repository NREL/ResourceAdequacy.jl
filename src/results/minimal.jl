struct MinimalResult <: ResultSpec end

# Single-period reliability results

struct SinglePeriodMinimalResult{
    N,T<:Period,P<:PowerUnit,E<:EnergyUnit,V<:AbstractFloat,
    SS<:SimulationSpec} <: SinglePeriodReliabilityResult{N,T,P,E,V,SS}

    lolp::LOLP{N,T,V}
    eue::EUE{E,N,T,V}
    simulationspec::SS

    function SinglePeriodMinimalResult{P}(
        lolp::LOLP{N,T,V},
        eue::EUE{E,N,T,V},
        simulationspec::SS) where {
            N,T,P<:PowerUnit,E,V,SS<:SimulationSpec}
        new{N,T,P,E,V,SS}(lolp, eue, simulationspec)
    end

end

LOLP(x::SinglePeriodMinimalResult) = x.lolp
EUE(x::SinglePeriodMinimalResult) = x.eue

mutable struct SinglePeriodNetworkMinimalResultAccumulator{
    N,T<:Period,P<:PowerUnit,E<:EnergyUnit,V<:Real}

    lol_count::Int
    eue::OnlineStats.Variance{OnlineStatsBase.EqualWeight}
    edgelabels::Vector{Tuple{Int,Int}}
    simulationspec::NonSequentialNetworkFlow

    function SinglePeriodNetworkMinimalResultAccumulator{}(
        simulationspec::NonSequentialNetworkFlow,
        sys::SystemStateDistribution{N,T,P,E,V}
    ) where {N,T,P,E,V}
        new{N,T,P,E,V}(0, Variance(), sys.interface_labels, simulationspec)
    end
end

function update!(acc::SinglePeriodNetworkMinimalResultAccumulator{N,T,P,E,Float64},
                 statematrix::Matrix{Float64},
                 flowmatrix::Matrix{Float64},
                 sink_idx::Int, n_regions::Int) where {N,T,P,E}

    if !all_load_served(statematrix, flowmatrix, sink_idx, n_regions)
        acc.lol_count += 1
        ns = NetworkState{N,T,P,E}(statematrix, flowmatrix, acc.edgelabels, n_regions)
        fit!(acc.eue, powertoenergy(droppedload(ns), N, T, P, E))
    else
        fit!(acc.eue, 0.)
    end

    return acc

end

function finalize(acc::SinglePeriodNetworkMinimalResultAccumulator{N,T,P,E,V}
                  ) where {N,T,P,E,V}

    nsamples = acc.simulationspec.nsamples

    lolp_μ = acc.lol_count/nsamples
    lolp_σ² = lolp_μ * (1-lolp_μ)

    eue_μ = acc.eue.μ
    eue_σ² = acc.eue.σ2

    return SinglePeriodMinimalResult{P}(
        LOLP{N,T}(lolp_μ, sqrt(lolp_σ²/nsamples)),
        EUE{E,N,T}(eue_μ, sqrt(eue_σ²/nsamples)),
        acc.simulationspec)

end

# Multi-period reliability results

struct MultiPeriodMinimalResult{
    N1,
    T1<:Period,
    N2,
    T2<:Period,
    P<:PowerUnit,
    E<:EnergyUnit,
    V<:AbstractFloat,
    ES<:ExtractionSpec,
    SS<:SimulationSpec
} <: MultiPeriodReliabilityResult{N1,T1,N2,T2,P,E,V,ES,SS}

    timestamps::StepRange{DateTime,T1}
    results::Vector{SinglePeriodMinimalResult{N1,T1,P,E,V,SS}}
    extractionspec::ES
    simulationspec::SS

    function MultiPeriodMinimalResult{}(
        timestamps::StepRange{DateTime,T1},
        results::Vector{SinglePeriodMinimalResult{N1,T1,P,E,V,SS}},
        extractionspec::ES,
        simulationspec::SS
    ) where {N1,T1,P,E,V,ES<:ExtractionSpec,SS<:SimulationSpec}
        n = length(timestamps)
        @assert n == length(results)
        @assert step(timestamps) == T1(N1)
        new{N1,T1,N1*n,T1,P,E,V,ES,SS}(
            timestamps, results, extractionspec, simulationspec)
    end
end

function MultiPeriodMinimalResult(
    dts::StepRange{DateTime,T}, results::Vector{R},
    extrspec::ES) where {T<:Period,R<:SinglePeriodMinimalResult, ES<:ExtractionSpec}

    simulationspec = results[1].simulationspec
    return MultiPeriodMinimalResult(dts, results, extrspec, simulationspec)

end

function MultiPeriodMinimalResult(
    sys::SystemModel{N1,T1,N2,T2,P,E,V},
    shortfalls::Matrix{V},
    extractionspec::ExtractionSpec,
    simulationspec::SimulationSpec
) where {N1, T1<:Period, N2, T2<:Period,
         P<:PowerUnit, E<:EnergyUnit, V<:Real}

    n_periods, n_samples = size(shortfalls)
    ps = zeros(V, n_periods)
    eues = zeros(V, n_periods)

    for i in 1:n_samples, t in 1:n_periods
        shortfall = shortfalls[t, i]
        if shortfall > 0
            ps[t] += 1.
            eues[t] += shortfall
        end
    end

    ps ./= n_samples
    eues ./= n_samples

    p_stderrs = sqrt.(ps.*(1.-ps) ./ n_samples)
    eue_stderrs = sqrt.(vec(var(shortfalls, 2, mean=eues))./n_samples)

    return MultiPeriodMinimalResult(
        sys.timestamps,
        SinglePeriodMinimalResult{P}.(
            LOLP{N1,T1}.(ps, p_stderrs),
            EUE{E,N1,T1}.(eues, eue_stderrs),
            simulationspec),
        extractionspec)

end

aggregator(::MinimalResult) = MultiPeriodMinimalResult

timestamps(x::MultiPeriodMinimalResult) = x.timestamps

function Base.getindex(x::MultiPeriodReliabilityResult,
                       dt::DateTime)
    idxs = searchsorted(x.timestamps, dt)
    if length(idxs) > 0
        return x.results[first(idxs)]
    else
        throw(BoundsError(x, dt))
    end
end

LOLE(x::MultiPeriodMinimalResult) = LOLE([LOLP(r) for r in x.results])
EUE(x::MultiPeriodMinimalResult) = EUE([EUE(r) for r in x.results])
