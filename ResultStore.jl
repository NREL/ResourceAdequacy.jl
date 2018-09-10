struct ResultStore
    firmcapacities::Vector{Float64}
    metrics::Vector{Float64}

    function ResultStore(fcs::Vector{Float64}, ms::Vector{Float64})
	@assert length(fcs) == length(ms)
        sortorder = sortperm(ms)
        permute!(fcs, sortorder)
        permute!(ms, sortorder)
        new(fcs, ms)
    end

end

function update!(x::ResultStore,
                 newfirmcapacities::Vector{Float64},
                 newmetrics::Vector{Float64})

    @assert length(newfirmcapacities) == length(newmetrics)

    append!(x.firmcapacities, newfirmcapacities)
    append!(x.metrics, newmetrics)

    neworder = sortperm(tuple.(x.metrics, .-x.firmcapacities))
    permute!(x.firmcapacities, neworder)
    permute!(x.metrics, neworder)

    foreach((fc, m) -> println(fc, " MW => ", m),
            x.firmcapacities, x.metrics)
            
    issorted(x.firmcapacities, rev=true) || error(
        "Metric monotonicity violation, " *
        "try increasing sample size to reduce sampling error")

    return x

end

Base.length(x::ResultStore) = length(x.firmcapacities)

writeresult(path::String, result::ResultStore) =
    writecsv(path, [result.firmcapacities result.metrics])
