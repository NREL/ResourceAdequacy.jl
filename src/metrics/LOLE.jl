# Loss-of-Load Expectation

type LOLE{N1,P1<:Period,N2,P2<:Period,V<:AbstractFloat} <: ReliabilityMetric{V}
    val::V
    stderr::V

    function LOLE{N1,P1,N2,P2}(val::V, stderr::V) where {N1,P1<:Period,N2,P2<:Period,V<:AbstractFloat}
        (val >= 0) || error("$val is not a valid occurence expectation")
        (stderr >= 0) || error("$stderr is not a valid standard error")
        new{N1,P1,N2,P2,V}(val, stderr)
    end

end

Base.show(io::IO, x::LOLE{N1,P1,N2,P2}) where {N1,P1,N2,P2} =
    print(io, "LOLE = ", val(x),
          stderr(x) > 0 ? "±"*string(stderr(x)) : "", " ",
          N1 == 1 ? "" : N1, unitsymbol(P1), "/",
          N2 == 1 ? "" : N2, unitsymbol(P2))

