# Expected Unserved Energy

struct EUE{E<:EnergyUnit,N,L,T<:Period,V<:Real} <: ReliabilityMetric{V}
    val::V
    stderr::V

    function EUE{E,N,L,T}(val::V, stderr::V) where {E<:EnergyUnit,N,L,T<:Period,V<:Real}
        (val >= 0) || error("$val is not a valid unserved energy expectation")
        (stderr >= 0) || error("$stderr is not a valid standard error")
        new{E,N,L,T,V}(val, stderr)
    end

end

Base.show(io::IO, x::EUE{E,N,L,T}) where {E,N,L,T} =
    print(io, "EUE = ", val(x),
          stderr(x) > 0 ? "±"*string(stderr(x)) : "", " ",
          unitsymbol(E), "/",
          N*L == 1 ? "" : N*L, unitsymbol(T))
