CapacityDistribution{T} = Distributions.Generic{T,Float64,Vector{T}}
CapacitySampler{T} = Distributions.GenericSampler{T, Vector{T}}

SumVariance{T} = OnlineStats.Series{
    Number,
    Tuple{OnlineStats.Sum{T},
          OnlineStats.Variance{OnlineStatsBase.EqualWeight}
}}

function mean_stderr(sv::SumVariance{T}) where T
    resultsum, _ = value(sv)
    return (resultsum, zero(T))
end

function mean_stderr(sv::SumVariance, nsamples::Int)
    samplesum, samplevar = value(sv)
    return (samplesum / nsamples, sqrt(samplevar / nsamples))
end

function searchsortedunique(a::AbstractVector{T}, i::T) where {T}
    idxs = searchsorted(a, i)
    length(idxs) == 0 && error("Element $i does not exist in $a")
    length(idxs) > 1 && error("Element $i occurs more than once in $a")
    return first(idxs)
end

function findfirstunique(a::AbstractVector{T}, i::T) where T
    i_idx = findfirst(a, i)
    i_idx > 0 || error("Element $i does not exist in $a")
    return i_idx
end

"""
Allocate each RNG on its own thread.
Note that the seed alone is not enough to enforce determinism: the number of
threads used will also affect results. For full reproducibility the thread
count should be constant between runs.
"""
function init_rngs(seed::UInt=rand(UInt))
    nthreads = Threads.nthreads()
    rngs = Vector{MersenneTwister}(nthreads)
    rngs_temp = randjump(MersenneTwister(seed),nthreads)
    Threads.@threads for i in 1:nthreads
        rngs[i] = copy(rngs_temp[i])
    end
    return rngs
end
