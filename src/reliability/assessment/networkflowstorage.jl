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

function all_load_served(A::Matrix{T}, B::Matrix{T}, sink::Int, n::Int) where T
    served = true
    i = 1
    while served && (i <= n)
        served = A[i, sink] == B[i, sink]
        i += 1
    end
    return served
end

function assess(params::NetworkFlowStorage, system::SystemDistribution{N,T,P,Float64}) where {N,T,P}

    systemsampler = SystemSampler(system)
    sink_idx = nv(systemsampler.graph) #Number of total nodes in the graph, sink is the last node
    source_idx = sink_idx-1 #Source is the second-to-last last node
    n = sink_idx-2 #Number of network nodes

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


    #In this section, use the generator parameters to determine the transition probabilities for each.
    generator_state_trans_matrix = Matrix{Float64}(length(system.gen_MTTR),2)
    for i in 1:length(system.gen_MTTR)
        #Create a state transition matrix from the MTTR and MTBF
        # 0-1, 1-0, the other transition probs are the complement of these values
        generator_state_trans_matrix[i,:] = [1./system.gen_MTTR[i], 1./system.gen_MTBF[i]]
    end

    initial_generator_ON_fraction = 0.8
    generator_state_vector = Int.(rand(length(system.gen_MTTR),1) .< initial_generator_ON_fraction*ones(length(system.gen_MTTR)) #Initialize a vector of ones (Generator ON) and zeros (Generator OFF)

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

        #Here, we need an update function to update each item based on previous timestep
        #Need some way to sum up all of the generation within an area. How should we do this?


        rand!(state_matrix, systemsampler)
        systemload, flow_matrix =
            LightGraphs.push_relabel!(flow_matrix, height, count, excess, active,
                          systemsampler.graph, source_idx, sink_idx, state_matrix)

        if !all_load_served(state_matrix, flow_matrix, sink_idx, n)

            # TODO: Save whether generator or transmission constraints are to blame?
            lol_count += 1
            #lol_sum += 0

            params.persist && push!(failure_states, FailureResult(state_matrix, flow_matrix, system.interface_labels, n))

        end

        #########################################################
        #Update generators (put this in a function??? Is there a more efficient way to do this??)
        for i in 1:length(system.gen_MTTR)
            if rand(1) .< generator_state_trans_matrix[i,generator_state_vector[i]+1]
                generator_state_vector[i] = Int.(~Bool(generator_state_vector[i]))
            else
                #Do nothing, keep state the same
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
