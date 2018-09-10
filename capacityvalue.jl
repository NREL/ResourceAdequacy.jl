function capacityvalue(metric::Symbol, targetmetrics::Vector{Float64}, tol::Float64,
                       maxcv::Float64, nparallel::Int, resultspath::String="")

    if metric == :lole
        metric_idx = 1
    elseif metric == :eue
        metric_idx = 2
    else
        error("Invalid metric $metric provided")
    end

    results = ResultStore(Float64[], Float64[])
    cvs = fill(NaN, length(targetmetrics))

    # Distribute initial solve points uniformly across the feasible range
    firmcapacities = collect(0:(maxcv/(nparallel-1)):maxcv)

    while any(isnan, cvs)
        metrics = getindex.(pmap(runsimulation, firmcapacities), metric_idx)
        update!(results, firmcapacities, metrics)
        prepnextiter!(cvs, firmcapacities, results, targetmetrics, tol)
    end

    resultspath != "" && writeresult(resultspath, results)

    return cvs

end

runsimulation(firm::Float64) =
    runsimulation(SimulationParams(0.1, 0.1, 0., 0.2, 10_000, "", firm))

"""
Store any firm capacities that correspond to a targetmetric as cvs, and 
allocate next round of firm capacities based on remaining unsolved target metrics
"""
function prepnextiter!(cvs::Vector{Float64}, firmcapacities::Vector{Float64},
                       results::ResultStore,
                       targetmetrics::Vector{Float64}, tol::Float64)

    n_results = length(results)
    n_targets = length(targetmetrics)
    result_idx = 2
    target_idx = 1

    nextranges = Set{Tuple{Float64,Float64}}()

    foreach((fc, m) -> println(fc, " MW => ", m), results.firmcapacities, results.metrics)
    while target_idx <= n_targets

        observedmetric = results.metrics[result_idx]
        targetmetric = targetmetrics[target_idx]

        if !isnan(cvs[target_idx])
            target_idx += 1
        elseif targetmetric < observedmetric - tol
            # observedmetric is the upper limit of the range containing target,
            # consider that range in the next iteration
            push!(nextranges, (results.firmcapacities[result_idx], results.firmcapacities[result_idx-1]))
            target_idx += 1
        elseif observedmetric - tol <= targetmetric <= observedmetric + tol
            # observed within tolerance of target, move to next target
            cvs[target_idx] = results.firmcapacities[result_idx]
            target_idx += 1
        else # observedmetric + tol  < targetmetric
            # observed is below remaining targets, move to next result observation
            result_idx += 1
        end

    end
    println([targetmetrics cvs])

    # Allocate n_results new tasks across nextranges
    n_ranges = length(nextranges)

    qtnt, rmdr = divrem(n_results, n_ranges)
    n_solvepoints = fill(qtnt, n_ranges)
    n_solvepoints[1:rmdr] += 1
    nextranges = sort!(collect(nextranges))
    fc_idx = 1
    println([nextranges n_solvepoints])

    for ((lower, upper), n) in zip(nextranges, n_solvepoints)
        spacing = (upper - lower) / (n+1)
        prev_fc = lower
        for i in fc_idx:(fc_idx+n-1)
            firmcapacities[i] = prev_fc+spacing
            prev_fc += spacing
        end
        fc_idx += n
    end

end


