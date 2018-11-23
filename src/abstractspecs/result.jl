# TODO: Metaprogrammed documentation for metrics

# Metrics defined over multiple timesteps
for T in [LOLP, LOLE, EUE] # LOLF would go here too

    # Metric over all timesteps and regions
    (::Type{T})(::R) where {R<:Result} = 
        error("$T(::$R) not yet defined")

    # Metric over all timesteps and specific region
    (::Type{T})(::R, ::AbstractString) where {R<:Result} = 
        error("$T(::$R, region::AbstractString) not defined: $R may not " *
              "support regional results")

    # Metric over specific timestep range and all regions
    (::Type{T})(::R, ::DateTime, ::DateTime) where {R<:Result} =
        error("$T(::$R, start::DateTime, end::DateTime) not defined: $R " *
              "may not support timestep sub-interval results")

    # Metric over specific region and specific timestep range
    (::Type{T})(::R, ::DateTime, ::DateTime, ::AbstractString) where {R<:Result} =
        error("$T(::$R, start::DateTime, end::DateTime, " *
              "region::AbstractString) not defined: $R may not support " *
              "regional timestep sub-interval results")

end

# Metrics defined over single timesteps
for T in [LOLP, EUE]

    # Metric at a specific timestep over all regions
    (::Type{T})(::R, ::DateTime) where {R<:Result} = 
        error("$T(::$R, period::DateTime) not yet defined: $R may not support" *
              "timestep-specific results")

   
    # Metric at a specific timestep and region
    (::Type{T})(::R, ::DateTime, ::AbstractString) where {R<:Result} = 
        error("$T(::$R, period::DateTime, region::AbstractString) " *
              "not yet defined: $R may not support regional " *
              "timestep-specific results")

end

"""

    accumulator(::ExtractionSpec, ::SimulationSpec, ::ResultSpec,
                ::SystemModel, seed::UInt)::ResultAccumulator

Returns a `ResultAccumulator` corresponding to the provided `ResultSpec`.
"""
accumulator(::ExtractionSpec, ::SimulationSpec, ::T,
            ::SystemModel, seed::UInt) where {T<:ResultSpec} = 
    error("An `accumulator` method has not been defined for ResultSpec $T")


# TODO: Finalize `update!` API

# Might have multiple possible type signatures here (copper plate + regional +
# everything?), with default conversions from more-complicated-but-general
# methods to simpler ones

# Note: `unservedenergyperiods` can't be aggregated across regions if the
#       results are an LOLP... Need to figure that out

"""

    update!(::ResultAccumulator{V}, unservedenergy::V,
            unservedenergyperiods::V, t::Int)

Records a simulation sample or result for a single timestep `t` in the
provided `ResultAccumulator`.

Implementation note: This function should be thread-safe as it will
generally be parallelized across many time periods and/or samples during
simulations. This is commonly achieved by storing results in a
thread-specific temporary storage location in the `ResultAccumulator` struct
and then merging results from all threads during `finalize`.
"""
update!(::R,
        unservedenergy::V,
        unservedenergyperiods::V, t::Int) where {V, R <: ResultAccumulator{V}} =
    error("update! has not yet been defined for ResultAccumulator $A")

"""

    update!(acc::ResultAccumulator{V}, unservedenergy::Vector{V},
            unservedenergyperiods::Vector{V}, t::Int)

Store region-specific results for a single time period.
The default behaviour is just to aggregate the regional results together and
call the system-wide method, but `ResultAccumulator`s can define their own
specialized methods to use the raw disaggregated data instead.
"""
update!(acc::ResultAccumulator{V},
        unservedenergy::Vector{V}, 
        unservedenergyperiods::Vector{V}, t::Int) where V =
    update!(acc, sum(unservedenergy), _, t)

"""

    finalize(::ExtractionSpec, ::SimulationSpec, ::ResultAccumulator)::Result

Returns a `Result` corresponding to the provided `ResultAccumulator`.
"""
finalize(::ExtractionSpec, ::SimulationSpec, ::A) where {A <: ResultAccumulator} =
    error("finalize not defined for ResultAccumulator $A")

