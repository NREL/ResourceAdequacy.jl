struct NonSequentialCopperplate <: SimulationSpec end

issequential(::NonSequentialCopperplate) = false
ismontecarlo(::NonSequentialCopperplate) = false
iscopperplate(::NonSequentialCopperplate) = true

function assess!(acc::ResultAccumulator,
                 simulationspec::NonSequentialCopperplate,
                 sys::SystemStateDistribution{N,T,P,E},
                 t::Int) where {N,T,P,E}

    # Collapse net load
    netloadsamples = vec(sum(sys.loadsamples, 1) .- sum(sys.vgsamples, 1))
    netload = to_distr(netloadsamples)

    # Collapse regions
    # (hopefully already done during extraction, this approach is very slow)
    supply = sys.region_maxdispatchabledistrs[1]
    for i in 2:length(sys.region_maxdispatchabledistrs)
        supply = add_dists(supply, sys.region_maxdispatchabledistrs[i])
    end

    lolp_val, eul_val = assess(supply, netload)
    eue_val = powertoenergy(eul_val, N, T, P, E)
    update!(acc, eue, lolp, t)

end

function to_distr(vs::Vector)
    p = 1/length(vs)
    cmap = countmap(vs)
    return Generic(collect(keys(cmap)),
                   [p * w for w in values(cmap)])
end
