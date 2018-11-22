abstract type ReliabilityMetric{V<:Real} end

include("LOLP.jl")
include("LOLE.jl")
include("EUE.jl")

# Common getter methods
for T in [LOLP, LOLE, EUE]

    @eval val(x::($T)) = x.val
    @eval stderr(x::($T)) = x.stderr

    @eval Base.isapprox(x::M, y::M) where {M<:($T)} =
        isapprox(x.val, y.val) &&
        isapprox(x.stderr, y.stderr)

end

# Note: Result-specific constructor methods are defined
#       in abstractspecs/results.jl and results/*.jl
