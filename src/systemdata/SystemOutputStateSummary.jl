struct SystemOutputStateSummary{L, T<:Period, E<:EnergyUnit, V<:Real}
    lolp_system::V
    lolp_regions::Vector{V}
    eue_regions::Vector{V}
end
