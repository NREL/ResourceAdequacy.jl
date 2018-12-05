struct Minimal <: ResultSpec end

struct TemporalResultAccumulator{V,S,ES,SS} <: ResultAccumulator{V,S,ES,SS}

end

struct TemporalResult{
    N, # Number of timesteps simulated
    L, # Length of each timestep
    T <: Period, # Units of timestep duration
    E <: EnergyUnit, # Units for energy results
    V <: Real, # Numerical type of value data
    ES <: ExtractionSpec,
    SS <: SimulationSpec
} <: Result{N,L,T,V,ES,SS}

end

LOLE(x::TemporalResult)
LOLP(x::TemporalResult, dt::DateTime)

EUE(x::TemporalResult)
EUE(x::TemporalResult, dt::DateTime)

function accumulator(extractionspec::ExtractionSpec,
                     simulationspec::SimulationSpec,
                     resultspec::Temporal, sys::SystemModel{N,L,T,P,E,V},
                     seed::UInt) where {N,L,T,P,E,V}

    return TemporalResultAccumulator()

end

function update!(acc::TemporalResultAccumulator{V},
                 sample::SystemOutputStateSample, t::Int, i::Int) where {V}
    return

end

function update!(acc::TemporalResultAccumulator,
                 result::SystemOutputStateSummary, t::Int)
    return

end

function finalize(acc::MinimalResultAccumulator{V,<:SystemModel{N,L,T,P,E,V}}
                  ) where {N,L,T,P,E,V}

    return TemporalResult()

end
