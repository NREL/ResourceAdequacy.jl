# Expected Unserved Energy

struct EUE{N,L,T<:Period,E<:EnergyUnit,V<:Real} <: ReliabilityMetric{V}
    val::V
    stderr::V

    function EUE{N,L,T,E}(val::V, stderr::V) where {N,L,T<:Period,E<:EnergyUnit,V<:Real}
        (val >= 0) || error("$val is not a valid unserved energy expectation")
        (stderr >= 0) || error("$stderr is not a valid standard error")
        new{N,L,T,E,V}(val, stderr)
    end

end

Base.show(io::IO, x::EUE{N,L,T,E}) where {N,L,T,E} =
    print(io, "EUE = ", val(x),
          stderr(x) > 0 ? "±"*string(stderr(x)) : "", " ",
          unitsymbol(E), "/",
          N*L == 1 ? "" : N*L, unitsymbol(T))
