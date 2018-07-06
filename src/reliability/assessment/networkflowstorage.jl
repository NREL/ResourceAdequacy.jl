struct NetworkFlowStorage <: ReliabilityAssessmentMethod
    timesteps::Int #Number of iterations
    persist::Bool   #True: Save detailed results to HDF5 file
                    #False: Save just LOLE

#Constructor function
    function NetworkFlowStorage(timesteps::Int, persist::Bool=false)
        @assert timesteps > 0
        new(timesteps, persist)
    end
end

function all_load_served(params::NetworkFlowStorage, A::Matrix{T}, B::Matrix{T}, sink::Int, n::Int) where T
    served = true
    unserved_load_data = zeros(0,2)
    i = 1
    while i <= n
        served = A[i, sink] == B[i, sink]
        if !served
            unserved_load_data = [unserved_load_data; [i A[i, sink] - B[i, sink]]] #Appends a new row each time with the node at which unserved load occurred, and the amount unserved.
        end
        i += 1
    end
    return unserved_load_data
end

function assess(params::NetworkFlowStorage, system::SystemDistribution{N,T,P,Float64}) where {N,T,P}

    systemsampler = SystemSampler(system)
    sink_idx = nv(systemsampler.graph) #Number of total nodes in the graph, sink is the last node
    source_idx = sink_idx-1 #Source is the second-to-last last node
    n = sink_idx-2 #Number of network nodes
    n_gens = size(system.gen_dists,1) #Number of system dispatchable generators, not necessarily equal to the number of nodes/areas
    n_storage = size(system.storage_params,1) #Number of energy storage devices

    #Create a state matrix and initialize counters
    #state_matrix contains line flows, including flows between source to load to sink
    state_matrix = zeros(sink_idx, sink_idx)
    lol_count = 0
    lol_sum = 0.
    failure_states = FailureResult{Float64}[]

    flow_matrix = Array{Float64}(sink_idx, sink_idx)
    height = Array{Int}(sink_idx)
    count = Array{Int}(2*sink_idx+1)
    excess = Array{Float64}(sink_idx)
    active = Array{Bool}(sink_idx)


    ###########################################################################
    #In this section, use the generator parameters to determine the transition probabilities for each.
    #Create a state transition matrix from the MTTR and FOR
    # [OFF to ON, ON to OFF], the other transition probs are the complement of these values
    #MTBF = MTTR/FOR - MTTR
    generator_state_trans_matrix = Matrix{Float64}(n_gens,2)
    generator_state_trans_matrix = [1./system.gen_dists[:,3] 1./(system.gen_dists[:,3]./system.gen_dists[:,4] - system.gen_dists[:,3])]

    initial_generator_ON_fraction = 0.9 #Percentage of generators that begin online
    generator_state_vector = Int.(rand(n_gens,1) .< initial_generator_ON_fraction*ones(n_gens) #Initialize a vector of ones (Generator ON) and zeros (Generator OFF)
    ###########################################################################

    ###########################################################################
    #Initialize energy storage
    #Add a final column to determine the max this storage can provide in the next timestep
    #MaxPowerAvail X timestep = MaxEnergy X SOC
    #Since timestep = 1 unit, MaxPowerAvail = MaxEnergy X SOC, so long as it is less than MaxPowerCap
    storage_energy_tracker = [system.storage_params min.(system.storage_params[:,2].*system.storage_params[:,3],system.storage_params[:,1])] #Last column is min(energy X SOC, MaxPowerCap)
    ###########################################################################


    #Main Loop
    for i in 1:params.timesteps

        #Iterate over all timesteps
        #This is where load flow is performed, using the push_relabel! function

#=
        A) Generation--Markov Chain with state transition matrixF
            Start with initial generator distribution for first time step. Then update using Markov Chain
        B) Load--Design for multiple possibilities:
            1) Preallocate load data (backcast?)
            2) Extraction "windowing" method
        C) Variable Gen--Same as Load
        D) Storage--If used, reduce SOC by amount used in an hour. If charged, increase SOC.
                If neither, model some trickle losses.
=#

        #Create state_matrix from systemsampler
        rand!(state_matrix, systemsampler)

        #Sum all of the generation in each node, changes each time iteration
        generation_row_vector = zeros(1,sink_idx)
        for i in 1:n
            temp_index_vector = find(system.gen_dists[:,2].==i) #temporary vector of the indices of the generators at node i
            generation_row_vector[i] = sum(system.gen_dists[temp_index_vector].*generator_state_vector[temp_index_vector]) #Multiply the generator capacities by the state vector, which will either multiply by zero or one, then sum. Even though system.gen_dists is a matrix, the indices will extract the values in the first column, as Julia is column-major
        end

        #Sum all of the storage in each node, changes each time iteration
        storage_row_vector = zeros(1,sink_idx)
        for i in 1:n
            temp_index_vector = find(storage_energy_tracker[:,4].==i) #temporary vector of the indices of the storage devices at node i
            temp_var = storage_energy_tracker[:,5] #Temporarily extract just the final column containing the storage_energy_tracker matrix, which contains just the max available power for this timestep for each storage device
            storage_row_vector[i] = sum(temp_var[temp_index_vector]) #temp_index_vector indexes into the "fifth column" of storage_energy_tracker which contains the power capacity available at this timestep, calculated from SOC and MaxEnergy
        end

        #TODO: Currently overwriting the values in the state_matrix, but should probably just directly save these in SystemSampler
        state_matrix[source_idx,:] = generation_row_vector


        systemload, flow_matrix =
            LightGraphs.push_relabel!(flow_matrix, height, count, excess, active,
                          systemsampler.graph, source_idx, sink_idx, state_matrix)

        #Create a state_matrix with another node for storage
        state_matrix_storage = state_matrix
        state_matrix_storage = [state_matrix_storage zeros(sink_idx,1)]
        state_matrix_storage = [state_matrix_storage; ]

        #Functions
        if !all_gen_used
            #Charge some batteries, if possible

        end

        unserved_load_data = all_load_served(params, state_matrix, flow_matrix, sink_idx, n)
        if size(unserved_load_data,1) > 0 #If one instance of unserved load exists

            #First, we check if storage can fix the problem
            #Add available storage to generation section in state_matrix and re-run push_relabel!
            #state_matrix[source_idx,:] = state_matrix[source_idx,:] + storage_energy_tracker[:,1]


            # TODO: Save whether generator or transmission constraints are to blame?
            lol_count += 1
            #lol_sum += 0

            params.persist && push!(failure_states, FailureResult(state_matrix, flow_matrix, system.interface_labels, n))

        end

        #########################################################
        #Update generators based on their state transition matrix
        #TODO: Consider if there is a more efficient way to do this and perhaps put in its own function
        for i in 1:n_gens
            temp_bitarray = rand(1) .< generator_state_trans_matrix[i,generator_state_vector[i]+1] #This expression creates a BitArray that can't be used in IF logic
            if temp_bitarray[1] #Extract the first value which is a boolean
                generator_state_vector[i] = Int.(~Bool(generator_state_vector[i]))
            else
                #Do nothing, keep state the same
        end
        #########################################################

        #########################################################
        #Update storage (put this in a function??? Is there a more efficient way to do this??)
        if #storage_used
            #increase or decrease SOC accordingly as well as MaxPowerAvail
            storage_energy_tracker
        else
            #Decrease due to losses
            storage_energy_tracker
        end
        #########################################################



    end

    μ = lol_count/params.timesteps
    σ² = μ * (1-μ)
    # eue_val, E = to_energy(lol_sum/params.timesteps, P, N, T)
    eue_val, E = to_energy(Inf, P, N, T)

    detailed_results = FailureResultSet(failure_states, system.interface_labels)

    return SinglePeriodReliabilityAssessmentResult(
        LOLP{N,T}(μ, sqrt(σ²/params.timesteps)),
        EUE{E,N,T}(eue_val, 0.),
        params.persist ? detailed_results : nothing
    )

end
