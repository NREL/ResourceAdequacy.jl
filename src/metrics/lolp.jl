# Loss-of-Load Probability

type LOLP{N,P<:Period,V<:AbstractFloat} <: ReliabilityMetric{V}
    val::V
    stderr::V

    function LOLP{N,P}(val::V, stderr::V) where {N,P<:Period,V<:AbstractFloat}
        (0 <= val <= 1) || error("$val is not a valid probability")
        (stderr >= 0) || error("$stderr is not a valid standard error")
        new{N,P,V}(val, stderr)
    end

end

Base.show(io::IO, x::LOLP{N,P}) where {N,P<:Period} =
    print(io, "LOLP = ", val(x),
          stderr(x) > 0 ? "±"*string(stderr(x)) : "",
          "/", N == 1 ? "" : N, unitsymbol(P))

# Abstract methods

# Result for all regions and times
LOLP(::T) where {T<:Result} = 
    error("LOLP(::$T) not yet defined")

# Result for a specific time
LOLP(::T, ::DateTime) where {T<:Result} = 
    error("LOLP(::$T, period::DateTime) not yet defined")

# Result for a specific region
LOLP(::T, ::AbstractString) where {T<:Result} = 
    error("LOLP(::$T, region::AbstractString) not yet defined")

# Result for a specfic region and time
LOLP(::T, ::DateTime, ::AbstractString) where {T<:Result} = 
    error("LOLP(::$T, period::DateTime, region::AbstractString) " *
          "not yet defined")
