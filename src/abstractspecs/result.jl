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


"""

    savetimestepsample!(::ResultAccumulator, _)

Records the results from a single timestep in a single Monte Carlo sample in
the provided `ResultAccumulator`.

Implementation note: This function should be thread-safe as it
will generally be parallelized across many samples and/or timesteps during
simulations. This is commonly achieved by storing the results in a
thread-specific temporary storage location in the `ResultAccumulator` struct
and then merging results from all threads during `finalize`.
"""
savetimestepsample!(::ResultAccumulator, _) #TODO: Finalize API(s)
# Might have multiple possible type signatures here (copper plate + regional +
# everything?), with default conversions from more-complicated-but-general
# methods to simpler ones

"""

    savetimestepresult!(::ResultAccumulator, _)

Records the final results of the simulation for a single timestep in the
provided `ResultAccumulator`.

Implementation note: This function should be thread-safe as it will
generally be parallelized across many time periods during
simulations. This is commonly achieved by storing cross-timestep results in a
thread-specific temporary storage location in the `ResultAccumulator` struct
and then merging results from all threads during `finalize`.
"""
savetimestepresult!(::ResultAccumulator, _) #TODO: Finalize API(s)
# Might have multiple possible type signatures here (copper plate + regional +
# everything?), with default conversions from more-complicated-but-general
# methods to simpler ones

"""

    finalize(::ExtractionSpec, ::SimulationSpec, ::ResultAccumulator)::Result

Returns a `Result` corresponding to the provided `ResultAccumulator`.
"""
finalize(::ExtractionSpec, ::SimulationSpec, ::A) where {A <: ResultAccumulator} =
    error("finalize not defined for ResultAccumulator $A")

