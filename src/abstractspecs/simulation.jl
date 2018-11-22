"""

   iscopperplate(::SimulationSpec)::Bool

Defines whether the simulation method associated with `SimulationSpec` will
ignore disaggregated region and transmission data. Computational speedups are
generally possible when this is true.
"""
iscopperplate(::S) where {S <: SimulationSpec} = 
    error("iscopperplate not yet defined for SimulationSpec $S")


"""

    assess!(::ResultAccumulator, ::ExtractionSpec,
            ::SimulationSpec{Sequential}, ::SystemModel)

Run a full simulation of `SystemModel` according to `ExtractionSpec` and
`SimulationSpec`, storing the results in `ResultAccumulator`.
"""
assess!(::ResultAccumulator, ::ExtractionSpec, ::S, ::SystemModel
) where {S <: SimulationSpec{Sequential}} =
    error("assess! not yet defined for SimulationSpec $S")

"""

    assess!(::ResultAccumulator, ::SimulationSpec{NonSequential},
            ::SystemStateDistribution, t::Int)

Solve a `SystemStateDistribution` at timestep `t` of the simulation as
specified by `SimulationSpec`, storing the results in `ResultAccumulator`.
"""
assess!(::ResultAccumulator, ::S, ::SystemStateDistribution, t::Int
) where {S <: SimulationSpec{NonSequential}} =
    error("assess! not yet defined for SimulationSpec $S")

function assess(extractionspec::ExtractionSpec,
                simulationspec::SimulationSpec{NonSequential},
                resultspec::ResultSpec,
                system::SystemModel,
                seed::UInt=rand(UInt))

    acc = accumulator(extractionspec, simulationspec, resultspec, system, seed)
    statedistrs = extract(extractionspec, system, iscopperplate(simulationspec))

    Threads.@threads for (t, statedistr) in enumerate(statedistrs)
        assess!(acc, simulationspec, statedistr, t)
    end

    return finalize(extractionspec, simulationspec, acc)

end

function assess(extractionspec::ExtractionSpec,
                simulationspec::SimulationSpec{Sequential},
                resultspec::ResultSpec,
                system::SystemModel,
                seed::UInt=rand(UInt))

    acc = accumulator(extractionspec, simulationspec, resultspec, system, seed)

    Threads.@threads for i in 1:simulationspec.nsamples
        assess!(acc, extractionspec, simulationspec, system)
    end

    return finalize(acc)

end
