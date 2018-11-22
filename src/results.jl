### Result

"""
An abstract parent type for simulation results. When defining a new type
`R where {R <: Result}`, you must define methods for the following functions
/ constructors:

 - `LOLP`
 - `LOLE`
 - `EUE`

Check the documentation for each function for required type signatures. #TODO

You must also define the following allied types and their associated methods:

 - `S where {S <: ResultSpec}`
 - `A where {A <: ResultAccumulator}`

"""
abstract type Result{
    N, # Number of timesteps simulated
    L, # Length of each simulation timestep
    T <: Period, # Units of each simulation timestep
    P <: PowerUnit, # Units for reported power values
    E <: EnergyUnit, # Units for reported energy values
    V <: Real, # Numeric type of value data
    ES <: ExtractionSpec, # Prob. distr. extraction method for input time series
    SS <: SimulationSpec, # Type of simulation that produced the result
} end

# TODO: Metaprogrammed documentation for metrics

# Metrics defined over multiple timesteps
for T in [LOLP, LOLE, EUE] # LOLF would go here too

    # Metric over all timesteps and regions
    @eval ($T)(::R) where {R<:Result} = 
        error("$T(::$R) not yet defined")

    # Metric over all timesteps and specific region
    @eval ($T)(::R, ::AbstractString) where {R<:Result} = 
        error("$T(::$R, region::AbstractString) not yet defined")

    # Metric over specific timestep range and all regions
    @eval ($T)(::R, ::DateTime, ::DateTime) where {R<:Result} =
        error("$T(::$R, start::DateTime, end::DateTime) not yet defined")

    # Metric over specific region and specific timestep range
    @eval ($T)(::R, ::DateTime, ::DateTime, ::AbstractString) where {R<:Result} =
        error("$T(::$R, start::DateTime, end::DateTime, region::AbstractString) " *
              "not yet defined")

end

# Metrics defined over single timesteps
for T in [LOLP, EUE]

    # Metric at a specific timestep over all regions
    @eval ($T)(::R, ::DateTime) where {R<:Result} = 
        error("$T(::$R, period::DateTime) not yet defined")

   
    # Metric at a specific timestep and region
    @eval ($T)(::R, ::DateTime, ::AbstractString) where {R<:Result} = 
        error("$T(::$R, period::DateTime, region::AbstractString) " *
              "not yet defined")

end


### ResultSpec

"""
An abstract parent type for specifying how results should be stored. When
defining a new type `S where {S <: ResultSpec}, you must define methods for
the following functions:

 - `accumulator`

Check the documentation for each function for required type signatures.
"""
abstract type ResultSpec end

"""

    accumulator(::ExtractionSpec, ::SimulationSpec, ::ResultSpec,
                ::SystemModel, seed::UInt)::ResultAccumulator

Returns a `ResultAccumulator` corresponding to the provided `ResultSpec`.
"""
accumulator(::ExtractionSpec, ::SimulationSpec, ::T,
            ::SystemModel, seed::UInt) where {T<:ResultSpec} = 
    error("An `accumulator` method has not been defined for ResultSpec $T")


### ResultAccumulator

"""
An abstract parent type for accumulating simulation results. When defining
a new type `A where {A <: ResultAccumulator}`, you must define methods for
the following functions:

 - `savetimestepsample!` - for Monte Carlo simulations
 - `savetimestepresult!` - for time-partitioned analytical solutions
 - `finalize`

Check the documentation for each function for required type signatures.
"""
abstract type ResultAccumulator{
    ES <: ExtractionSpec,
    SS <: SimulationSpec
} end

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


# Concrete instantiations
include("results/minimal.jl")
include("results/network.jl")
