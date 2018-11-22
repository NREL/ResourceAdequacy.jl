# Expected Unserved Energy

type EUE{E<:EnergyUnit,N,P<:Period,V<:AbstractFloat} <: ReliabilityMetric{V}
    val::V
    stderr::V

    function EUE{E,N,P}(val::V, stderr::V) where {E<:EnergyUnit,N,P<:Period,V<:AbstractFloat}
        (val >= 0) || error("$val is not a valid unserved energy expectation")
        (stderr >= 0) || error("$stderr is not a valid standard error")
        new{E,N,P,V}(val, stderr)
    end

end

Base.show(io::IO, x::EUE{E,N,P}) where {E<:EnergyUnit,N,P<:Period} =
    print(io, "EUE = ", val(x),
          stderr(x) > 0 ? "±"*string(stderr(x)) : "", " ",
          unitsymbol(E), "/",
          N == 1 ? "" : N, unitsymbol(P))
