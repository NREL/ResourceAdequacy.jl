struct SequentialCopperplate <: SimulationSpec{Sequential}
    nsamples::Int

    function SequentialCopperplate(nsamples::Int)
        @assert nsamples > 0
        new(nsamples)
    end
end

iscopperplate(::SequentialCopperplate) = true

function assess_singlesequence!(
    shortfalls::AbstractVector{V},
    rng::MersenneTwister,
    extractionspec::Backcast, #TODO: Generalize
    simulationspec::SequentialCopperplate,
    resultspec::MinimalResult, #TODO: Generalize
    sys::SystemModel{N1,T1,N2,T2,P,E,V}
) where {N1,T1,N2,T2,P,E,V}

    # Initialize generator and storage state vector
    # based on long-run probabilities from period 1
    # Note: Could pre-allocate these once per thread?
    gens_available = Bool[rand(rng) < gen.μ /(gen.λ + gen.μ)
                         for gen in view(sys.generators, :, 1)]
    stors_available = Bool[rand(rng) < stor.μ / (stor.λ + stor.μ)
                          for stor in view(sys.storages, :, 1)]
    stors_energy = zeros(V, size(sys.storages, 1))
    DR_energy = zeros(V, size(sys.DR,1)) #TODO: Add DR to SystemModel, might not need to be the same size
    DR_events_tracker = zeros(V,0,2) #Type V, begins with no events, each event has two columns: [Energy, Payback Periods Remaining]

    # Main simulation loop
    for (t, (gen_set, stor_set)) in enumerate(zip(
        sys.timestamps_generatorset, sys.timestamps_storageset))

        netload = sum(view(sys.load, :, t)) - sum(view(sys.vg, :, t))

        dispatchable_gen_available = available_capacity!(
            rng, gens_available, view(sys.generators, :, gen_set))

        # println("Net Load: ", netload, "\t",
        #         "Dispatchable Gen: ", dispatchable_gen_available)
        residual_generation = dispatchable_gen_available - netload

        if residual_generation >= 0

            #Repay shifted DR first
            residual_generation = repay_DR!(residual_generation,...)

            # Charge to consume residual_generation
            charge_storage!(rng, stors_available, stors_energy,
                            residual_generation,
                            view(sys.storages, :, stor_set))

        else

            # Discharge to meet residual_generation shortfall
            shortfall = discharge_storage!(
                rng, stors_available, stors_energy,
                -residual_generation, view(sys.storages, :, stor_set))

            if shortfall > 0
            # Inject DR as a last resort
            inject_DR!(shortfall, sys.DR) #will have to change the "sys.DR"
            else
                #Skip over the function call
            end

            # Report remaining shortfall, if any
            shortfall > 0 && (shortfalls[t] = shortfall)

        end

    end
    function update_DR_tracking() #This is just a placeholder for now, add shortfalls if time limit expires
end

function available_capacity!(rng::MersenneTwister,
                             gen_availability::Vector{Bool},
                             generators::AbstractVector{DispatchableGeneratorSpec{T}}
                             ) where {T <: Real}

    capacity = zero(T)

    @inbounds for i in 1:length(gen_availability)
        gen = generators[i]
        if gen_availability[i]
            if rand(rng) > gen.λ # Unit doesn't fail, count capacity
                capacity += gen.capacity
            else # Unit fails, ignore its capacity
                gen_availability[i] = false
            end
        else
            if rand(rng) < gen.μ # Unit is repaired, count its capacity
                gen_availability[i] = true
                capacity += gen.capacity
            end
        end
    end

    return capacity

end

function charge_storage!(rng::MersenneTwister,
                         stors_available::Vector{Bool},
                         stors_energy::Vector{T},
                         surplus::T,
                         stors::AbstractVector{StorageDeviceSpec{T}}
                         ) where {T <: Real}

    # TODO: Replace with strategic charging
    # TODO: Stop assuming hourly periods

    for (i, stor) in enumerate(stors)

        stors_energy[i] *= stor.decayrate

        if stors_available[i]

            if rand(rng) > stor.λ # Unit doesn't fail

                max_charge = stor.energy - stors_energy[i]

                if surplus > max_charge # Fully charge

                    stors_energy[i] = stor.energy
                    surplus -= max_charge

                else # Partially charge

                    stors_energy[i] += surplus
                    return

                end

            else # Unit fails
                stors_available[i] = false
            end

        else

            if rand(rng) < stor.μ # Unit is repaired
                stors_available[i] = true

                max_charge = stor.energy - stors_energy[i]

                if surplus > max_charge # Fully charge

                    stors_energy[i] = stor.energy
                    surplus -= max_charge

                else # Partially charge

                    stors_energy[i] += surplus
                    return

                end

            end
        end
    end

end

function discharge_storage!(rng::MersenneTwister,
                            stors_available::Vector{Bool},
                            stors_energy::Vector{T},
                            shortfall::T,
                            stors::AbstractVector{StorageDeviceSpec{T}}
                            ) where {T <: Real}

    # TODO: Replace with strategic discharging
    # TODO: Stop assuming hourly periods

    for (i, stor) in enumerate(stors)

        stors_energy[i] *= stor.decayrate

        if stors_available[i]

            if rand(rng) > stor.λ # Unit doesn't fail

                max_discharge = stors_energy[i]

                if shortfall > max_discharge # Fully discharge

                    stors_energy[i] = 0
                    shortfall -= max_discharge

                else # Partially discharge

                    stors_energy[i] -= shortfall
                    return zero(T)

                end

            else # Unit fails
                stors_available[i] = false
            end

        else

            if rand(rng) < stor.μ # Unit is repaired

                stors_available[i] = true

                max_discharge = stors_energy[i]

                if shortfall > max_discharge # Fully discharge

                    stors_energy[i] = 0
                    shortfall -= max_discharge

                else # Partially discharge

                    stors_energy[i] -= shortfall
                    return zero(T)

                end

            end
        end
    end

    return shortfall

end

function repay_DR!(surplus::T, DR_events_tracker::Array{U,2}  ) where {T <: Real, U <: Real}

    #TODO: THESE INPUTS ABOVE ARE JUST PLACEHOLDERS

    #Sort by periods remaining to repay
    #TODO: Should we sort by energy amount also?
    DR_events_tracker = sortslices(DR_events_tracker, dims=1,lt=(x,y)->isless(x[2],y[2]))

    if surplus > 0
        for (i, event) in enumerate(DR_events_tracker[:,2])
            if surplus >= DR_events_tracker[i,1]
                surplus -= DR_events_tracker[i,1]
                DR_events_tracker[i,:] = [0 0] #Event has been repayed, timer reset
            else #surplus < DR_events_tracker[i,1]
                DR_events_tracker[i,1] -= surplus
                surplus = 0
                return(zero(T))
            end
        end
    else #surplus <= 0
        #Do Nothing
        #Shouldn't be possible as this function is called only if surplus > 0
    end
end

 function inject_DR!(energy_required, available_DR
                          ) where {T <: Real}

    #Sort by payback time
    available_DR = sortslices(available_DR, dims=1,lt=(x,y)->isless(x[2],y[2]))
    DR_tracking = 0
    for (i, DR_energy) in enumerate(available_DR)
        if energy_required >= DR_energy
            DR_tracking += DR_energy
            energy_required -= DR_energy
            available_DR[i,1] = 0 #Might have no reason to keep track of this
        else #energy_required <= available_DR[]
            DR_tracking += energy_required
            available_DR[i,1] -= energy_required #Might not need to track this
            energy_required = 0

            #End the loop here
            return DR_tracking, energy_required
        end
    end

    return DR_tracking, energy_required
  end
