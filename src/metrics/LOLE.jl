# Loss-of-Load Expectation

struct LOLE{N,L,T<:Period,V<:Real} <: ReliabilityMetric{V}
    val::V
    stderr::V

    function LOLE{N,L,T}(val::V, stderr::V) where {N,L,T<:Period,V<:Real}
        (val >= 0) || error("$val is not a valid occurence expectation")
        (stderr >= 0) || error("$stderr is not a valid standard error")
        new{N,L,T,V}(val, stderr)
    end

end

Base.show(io::IO, x::LOLE{N,L,T}) where {N,L,T} =
    print(io, "LOLE = ", val(x),
          stderr(x) > 0 ? "±"*string(stderr(x)) : "", " ",
          L == 1 ? "" : L, unitsymbol(T), "/",
          N*L == 1 ? "" : N*L, unitsymbol(T))

