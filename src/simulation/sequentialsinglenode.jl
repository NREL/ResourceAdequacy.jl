using MathProgBase
using Clp
using RowEchelon
using DataFrames
using CSV
using JLD
#using ResourceAdequacy

#Inputs:
#Not running the assess function, just running the script
#Sequential one-node system (ignoring transmission)

#region_labels = ["A","B","C"]

# gen_dists = [Generic([2., 3], [.4, .6]),
#              Generic([2., 3], [.4, .6]),
#              Generic([2., 3], [.4, .6])]

#Generators: [MaxCap Node/Area MTTR FOR]
gen_distributions_sequential = [100 1 2 0.05;
             200 1 3 0.1;
             300 1 3 0.1;
             1000 1 2 0.08;
             800 1 3 0.05]

###############################################################################
#Data manipulation
temp_gen = CSV.read("C:/Users/aklem/Documents/GitHub/ResourceAdequacy/gen.csv",rows_for_type_detect = 200) #Default is rows_for_type_detect = 100 which causes an error
gen_distributions_sequential = Array(temp_gen[1:96,[11,1,28,26]]) #Extract rated power,Extracting arbitrary column to replace with node/area, MTTR (Hrs), FOR
gen_distributions_sequential[:,2] = ones(size(gen_distributions_sequential,1)) #All in node one for this work
#Need to scale the rated power accordingly

#Import the load data
load = CSV.read("C:/Users/aklem/Documents/GitHub/ResourceAdequacy/LA_residential_enduses.csv")
load_matrix = collect(Missings.replace(Array(load[:,2:9]),0))
heating_loads = load_matrix[:,3]
cooling_loads = load_matrix[:,4]
heating_load_year_max = maximum(heating_loads)
cooling_load_year_max = maximum(cooling_loads)
total_load = sum(load_matrix,2)

heating_loads_repay_time = 4
cooling_loads_repay_time = 2
usable_DR_fraction = 0 #Participation factor [0,1]

#Initialize timesteps and MC iterations
MonteCarloIterations = 1
timesteps = size(load_matrix,1)

#Add in solar and wind
solar_power = CSV.read("C:/Users/aklem/Documents/GitHub/ResourceAdequacy/DAY_AHEAD_pv.csv",rows_for_type_detect = 200)
wind_power = CSV.read("C:/Users/aklem/Documents/GitHub/ResourceAdequacy/DAY_AHEAD_wind.csv",rows_for_type_detect = 200)
solar_power = collect((Missings.replace(Array(solar_power[1:end,5:29]),0)))
wind_power = collect(Missings.replace(Array(wind_power[:,5:8]),0))

#These will be used in the main loop
solar_sample = zeros(size(solar_power,2))
wind_sample = zeros(size(wind_power,2))

#Determine wind and solar max values for each farm in the data, and their capacity factors
solar_maxima = findmax(solar_power,1)[1]
wind_maxima = findmax(wind_power,1)[1]

#Timesteps can be interchanged here because the load, solar, and wind are all yearly data with a "leap" day
solar_cap_factor = sum(solar_power,1)./timesteps./12./solar_maxima
wind_cap_factor = sum(wind_power,1)./timesteps./12./wind_maxima

#In order to properly scale the sizes of the capacities
dispatchable_gen_capacity = sum(gen_distributions_sequential[:,1].*(1-0*gen_distributions_sequential[:,4]))
solar_capacity = sum(solar_maxima.*solar_cap_factor)
wind_capacity = sum(wind_maxima.*wind_cap_factor)

desired_solar_fraction = 0.0
desired_wind_fraction = 0.0
gen_margin = 0.0

gen_scale_factor = (1 + gen_margin)*(1 - desired_solar_fraction - desired_wind_fraction)*maximum(total_load)./dispatchable_gen_capacity
solar_scale_factor = desired_solar_fraction*(1 + gen_margin)*maximum(total_load)./solar_capacity
wind_scale_factor = desired_wind_fraction*(1 + gen_margin)*maximum(total_load)./wind_capacity

#Overwrite previous generation values with the adjusted scalings
gen_distributions_sequential[:,1] *= gen_scale_factor
solar_power *= solar_scale_factor
wind_power *= wind_scale_factor

#End of load manipulation
###############################################################################
#Storage Parameters (Max Power, Max Energy Cap, Initial SOC, Node/Area)
storage_params = [1 4 0.9 1;
                 2 0.5 0.9 1;
                 0.5 0.5 0.9 1
                 1 2 0.9 1
                 2 1 0.9 1]
storage_params = [0 0 0 1]

#Demand Response Parameters (Power, Shiftable Periods, Time periods until payback is required, Node/Area)
# DR_params = [1 2 8 1
#             2 4 Inf 1]

#TODO: Currently hardcoded in as one T unit (e.g 1 hour, 1 minute). Fix this.

    n = 1 #One node
    n_gens = size(gen_distributions_sequential,1) #Number of system dispatchable generators, not necessarily equal to the number of nodes/areas
    n_storage = size(storage_params,1) #Number of energy storage devices

    lol_count = 0
    lol_sum = 0.
    lol_count_with_storage = 0

#Intialize outputs
output_data = zeros(timesteps,20)
UnservedLoadData = zeros(timesteps,MonteCarloIterations)
AvailableGenCap = zeros(timesteps,MonteCarloIterations)
DR_Injected = zeros(timesteps,MonteCarloIterations)
RunTime = zeros(MonteCarloIterations)

###########################################################################
#In this section, use the generator parameters to determine the transition probabilities for each.
#Create a state transition matrix from the MTTR and FOR
# [OFF to ON, ON to OFF], the other transition probs are the complement of these values
#MTBF = MTTR/FOR - MTTR
generator_state_trans_matrix = Matrix{Float64}(n_gens,2)
generator_state_trans_matrix = [1./gen_distributions_sequential[:,3] 1./(gen_distributions_sequential[:,3]./gen_distributions_sequential[:,4] - gen_distributions_sequential[:,3])]

initial_generator_ON_prob_vector = zeros(size(generator_state_trans_matrix,1))
for i in 1:size(generator_state_trans_matrix,1)
    temp = rref([-generator_state_trans_matrix[i,1] generator_state_trans_matrix[i,2] 0; 1 1 1]) #Will return a matrix [1 0 ProbOFF; 0 1 ProbON]
    initial_generator_ON_prob_vector[i] = temp[2,3]
    if isnan(temp[2,3])
        initial_generator_ON_prob_vector[i] = 1
    else
    end
end

generator_state_vector = Int.(rand(n_gens,1) .< initial_generator_ON_prob_vector.*ones(n_gens,1)) #Initialize a vector of ones (Generator ON) and zeros (Generator OFF)
###########################################################################

###########################################################################
#Initialize energy storage
#Add two columns to the input data to determine the max power this storage can provide in the next timestep and the max it can receive
#MaxPowerDischarge X timestep = MaxEnergy X SOC
#MaxPowerCharge X timestep = MaxEnergy X (1-SOC)
#Since timestep = 1 unit, MaxPowerAvail = MaxEnergy X SOC, so long as it is less than MaxPowerCap. Thus, we take the min
storage_energy_tracker = [
    storage_params min.(
        storage_params[:,2] .* storage_params[:,3], storage_params[:,1]
    ) min.(
        storage_params[:,2] .* (1 - storage_params[:,3]), storage_params[:,1]
    )] #Last column is min(energy X SOC, MaxPowerCap)
storage_energy_tracker = [storage_energy_tracker zeros(n_storage,1)] #The last column is to keep track of updates that will be implemented at the end of the loop. Originally, a vector was used, but with the sorting process that occurs, this can get difficult to manage
#storage_energy_tracker = [MaxPower MaxEnergy SOC Node AvailDischargeCap AvailChargeCap PowerUsedThisTimestep]
#TODO: Make storage_energy_tracker into a structure? Probably easier than a matrix
###########################################################################

###########################################################################
#Initialize DR here
#Create an empty array that will be filled in as necessary during each timestep
DR_energy_tracker = zeros(0,2) #[Max Power Cap (Depends on minimum of Heating or Cooling maxes and the available energy), Remaining Periods to pay back] Not keeping track of energy as  it is basically the same as power


###########################################################################


#TODO: Put this into a structure?
#Create the first n rows
gen_matrix = eye(n)
storage_discharge_matrix = eye(n)
DR_inject_matrix = eye(n)
unserved_load_matrix = eye(n)
storage_charge_matrix = eye(n)
DR_payback_matrix = eye(n)



#Generate the A matrix, rather the matrix of coefficients for linear optimization
#Group the elements together in the matrix i.e. [g1;g2,...;s1;s2;...;d1;d2;...]
A = [gen_matrix storage_discharge_matrix DR_inject_matrix unserved_load_matrix -storage_charge_matrix -DR_payback_matrix]
#Add DR constraint (DR injected plus DR repayed needs to be less than max DR power available)
A = [A; [0 0 1 0 0 1]]


#Create the objective vector. Initialize the length as the width of A (the number of objective variables)
gen_cost = 0
storage_discharge_cost = 1
DR_inject_cost = 2
unserved_load_cost = 1000
storage_charge_cost = -0.5
DR_payback_cost = -1.5

c = [gen_cost*ones(n); storage_discharge_cost*ones(n);
DR_inject_cost*ones(n); unserved_load_cost*ones(n); storage_charge_cost*ones(n);
DR_payback_cost*ones(n);
]

gen_lower_limits = zeros(n)
storage_discharge_lower_limits = zeros(n)
DR_inject_lower_limits = zeros(n)
unserved_load_lower_limits = zeros(n)
storage_charge_lower_limits = zeros(n)
DR_payback_lower_limits = zeros(n) #This may change during simulation as their "time limits" expire

for MCI in 1:MonteCarloIterations
    tic()
    for i in 1:timesteps
    #for i in 1:100

    ###############################################################################
        #Uncomment this if using the 5-min renewable generation data
        # #TODO: Make this a method in the SystemSampler function, perhaps?
        # #Sample solar and wind generation values
        # for j in 1:length(solar_sample)
        #     solar_sample[j] = solar_power[(i-1)*12 + rand(1:12),j]
        # end
        # for j in 1:length(wind_sample)
        #     wind_sample[j] = wind_power[(i-1)*12 + rand(1:12),j]
        # end

        #Uncomment this if using the hourly renewable generation data
        # for j in 1:length(solar_sample)
        #     solar_sample[j] = solar_power[i,j]
        #     wind_sample[j] = wind_power[i,j]
        # end
        solar_sample = solar_power[i,:]
        wind_sample = wind_power[i,:]
    ###############################################################################

        #For our case, lower bound = upper bound, and Ax = Load Demanded, not zero
        #Thus, Ax = lb = ub = Load
        lb = [total_load[i]; 0]
        ub = [total_load[i]; heating_load_year_max + cooling_load_year_max]
        storage_energy_tracker[:,7] = zeros(n_storage,1) #Resets the vector of the amount of power used by each storage device. It will be used to update their SOCs at the end of the outermost loop. Needs to be re-"zeroed" at the start of each loop

        temp_vector = zeros(size(DR_energy_tracker,1),1)
        for j in 1:size(DR_energy_tracker,1)
            if DR_energy_tracker[j,2] != 0
                temp_vector[j] = 0
            else #if it = 0, thus hitting its time limit
                temp_vector[j] = DR_energy_tracker[j,1]
            end
        end
        DR_payback_lower_limits = sum(temp_vector) #TODO: This would be different if it wasn't a one-node system

        #Sum all of the generation in each node, changes each time iteration
        gen_upper_limits = zeros(n)
        for j in 1:n
            temp_index_vector = find(gen_distributions_sequential[:,2].==j) #temporary vector of the indices of the generators at node j
            gen_upper_limits[j] = sum(gen_distributions_sequential[temp_index_vector].*generator_state_vector[temp_index_vector]) + sum(solar_sample) + sum(wind_sample) #Multiply the generator capacities by the state vector, which will either multiply by zero or one, then sum. Then add solar and wind Even though gen_dists is a matrix, the indices will extract the values in the first column, as Julia is column-major
        end

        #Sum all of the storage in each node, changes each time iteration
        storage_discharge_upper_limits = zeros(n)
        for j in 1:n
            temp_index_vector = find(storage_energy_tracker[:,4].==j) #temporary vector of the indices of the storage devices at node j
            temp_var = storage_energy_tracker[:,5] #Temporarily extract just the fifth column containing the storage_energy_tracker matrix, which contains just the max available power for this timestep for each storage device, as the max avail power may change each iteration
            storage_discharge_upper_limits[j] = sum(temp_var[temp_index_vector]) #temp_index_vector indexes into the "fifth column" of storage_energy_tracker which contains the power capacity available at this timestep, calculated from SOC and MaxEnergy
        end

        #TODO: Should we change this? Right now, it is not affected by the maximum load limitations
        #TODO: Ensure that the amount that can be injected doesn't interfere with the payback of DR
        DR_inject_upper_limits = usable_DR_fraction*(heating_loads[i] + cooling_loads[i])

        unserved_load_upper_limits = total_load[i] #Can't be any higher than the laod

        storage_charge_upper_limits = zeros(n)
        for j in 1:n
            temp_index_vector = find(storage_energy_tracker[:,4].==j) #temporary vector of the indices of the storage devices at node j
            temp_var = storage_energy_tracker[:,6] #Temporarily extract just the sixth column containing the storage_energy_tracker matrix, which contains just the max available power for this timestep for each storage device, as the max avail power may change each iteration
            storage_charge_upper_limits[j] = sum(temp_var[temp_index_vector]) #temp_index_vector indexes into the "sixth column" of storage_energy_tracker which contains the power capacity available at this timestep, calculated from SOC and MaxEnergy
        end

        DR_payback_upper_limits = sum(DR_energy_tracker[:,1])

        # gen_upper_limits
        # storage_discharge_upper_limits
        # DR_inject_upper_limits
        # unserved_load_upper_limits
        # storage_charge_upper_limits
        # DR_payback_upper_limits
        # #transmission_upper_limits, determined above

        l = [gen_lower_limits; storage_discharge_lower_limits; DR_inject_lower_limits;
            unserved_load_lower_limits; storage_charge_lower_limits;
            DR_payback_lower_limits]

        u = [gen_upper_limits; storage_discharge_upper_limits; DR_inject_upper_limits;
            unserved_load_upper_limits; storage_charge_upper_limits;
            DR_payback_upper_limits]


        solver = ClpSolver()
        solution = linprog(c,A,lb,ub,l,u,solver)


        #Make updates based on optimization outputs
        solution_vector = solution.sol

        gen_used = solution_vector[1:n]
        storage_discharged = solution_vector[n+1:2*n]
        storage_charged = solution_vector[4*n+1:5*n]
        unserved_load = solution_vector[3*n+1:4*n]
        DR_injected = solution_vector[2*n+1:3*n]
        DR_repayed = solution_vector[5*n+1:6*n]

        #TODO: This will be updated when more than one node is added as storage can be used at one node and charged at another
        @assert (storage_discharged==zeros(n)) | (storage_charged==zeros(n)) #They can't both be used since we have a one-node system

        #@assert (DR_injected==zeros(n)) | (DR_repayed==zeros(n))

        if storage_discharged != zeros(size(storage_discharged,1))

            sortrows(storage_energy_tracker, by=x->(x[4],x[5])) #Sort it by node, then by max discharge power available
            temp_vector = storage_energy_tracker[:,7] #Temporaily extract the "update" column from storage_energy_tracker. It will be merged back in.

            for j in 1:n #loop through all the nodes (not the sink or source)
                Reqd_power_node_j = storage_discharged[j]
                temp_indices = find(storage_energy_tracker[:,4].==j) #Find indices of node j

                #Go through the storage devices at node j, using the lowest-power one first (hopefully, we have a proof for this)
                for k in 1:length(temp_indices)
                    #Note that it should be impossible for Reqd_power_node_j > sum(storage_energy_tracker[temp_indices[k],5]), this loop should pose no problems

                    if Reqd_power_node_j > storage_energy_tracker[temp_indices[k],5]
                        #Store the max power this device can provide, then adjust the power required
                        temp_vector[temp_indices[k]] = storage_energy_tracker[temp_indices[k],5]
                        Reqd_power_node_j = Reqd_power_node_j - storage_energy_tracker[temp_indices[k],5]

                    else
                        temp_vector[temp_indices[k]] = Reqd_power_node_j
                        break #exit this FOR loop as we have already determined which storage devices will provide the requirement
                    end
                end
            end
            storage_energy_tracker[:,7] = temp_vector #Merge the changes back in

        elseif storage_charged != zeros(size(storage_charged,1))

            #Charge storage in tracker
            sortrows(storage_energy_tracker, by=x->(x[4],x[6])) #Sort it by node, then by max charge power available
            temp_vector = storage_energy_tracker[:,7]
            for j in 1:n
                Served_power_node_j = storage_charged[j] #The "load" provided by the storage devices
                temp_indices = intersect(find(storage_energy_tracker[:,4].==j),find(storage_energy_tracker[:,7].==0)) #Find indices of node j that were also not discharged in earlier steps
                #Go through the storage devices at node j, use lowest-power one first
                for k in 1:length(temp_indices)

                    #It should be impossible for Served_power_node_j > sum(storage_energy_tracker[temp_indices[k],6]), shouldn't pose a problem
                    if Served_power_node_j > storage_energy_tracker[temp_indices[k],6]

                        #Store the max power this device can provide (will be updated later), then adjust the remaining power that can be served
                        temp_vector[temp_indices[k]] = -storage_energy_tracker[temp_indices[k],6] #Negative means it will be charged
                        Served_power_node_j = Served_power_node_j - storage_energy_tracker[temp_indices[k],6]
                    else
                        temp_vector[temp_indices[k]] = -Served_power_node_j #Negative means it will be charged
                        break #exit the FOR loop as we have already determined which storage devices will provide the requirement
                    end
                end
            end
            storage_energy_tracker[:,7] = temp_vector #Merge the changes back in

        else
            #Do nothing. We will account for losses outside of the IF statement as leakage losses should occur no matter what
        end


        #TODO: This will have to be updated when more nodes are added as DR can be injected at one node, and received at another
        sortrows(DR_energy_tracker, by=x->(x[2]))
        if DR_injected != zeros(size(DR_injected,1))
            #Add a row
            if DR_injected[1] <= usable_DR_fraction*heating_loads[i]
                DR_energy_tracker = [DR_energy_tracker; [DR_injected heating_loads_repay_time]]
            else #DR <= usable_DR_fraction*(heating_loads + cooling_loads)
                DR_energy_tracker = [DR_energy_tracker; [usable_DR_fraction*heating_loads[i] heating_loads_repay_time]; [(DR_injected - usable_DR_fraction*heating_loads[i]) cooling_loads_repay_time]]
            end

        elseif DR_repayed != zeros(size(DR_injected,1))
            #TODO: Make this go over all nodes in the future
            temp = DR_repayed[1]
            for j in 1:size(DR_energy_tracker,1)
                if temp > DR_energy_tracker[j,1]
                    temp = temp - DR_energy_tracker[j,1]
                    DR_energy_tracker[j,:] = [0 0] #Pay back this item, move on to the next
                else
                    DR_energy_tracker[j,1] = DR_energy_tracker[j,1] - temp
                    break
                end
            end
        else
            #Do nothing
        end


        #########################################################
        #Update generators based on their state transition matrix
        #TODO: Consider if there is a more efficient way to do this and perhaps put in its own function
        for j in 1:n_gens
            temp_bitarray = rand(1) .< generator_state_trans_matrix[j,generator_state_vector[j]+1] #This expression creates a BitArray that can't be used in IF logic
            if temp_bitarray[1] #Extract the first value which is a boolean
                generator_state_vector[j] = Int.(~Bool(generator_state_vector[j]))
            else
                #Do nothing, keep state the same
            end
        end
        #########################################################

        #########################################################
        #Update storage SOC (put this in a function??? Is there a more efficient way to do this??)
        #TODO: Should we account for losses during charge/discharge? Right now it's 100% round-trip efficiency
        for j in 1:n_storage
            if storage_energy_tracker[j,7] > 0
                #Discharge
                storage_energy_tracker[j,3] = storage_energy_tracker[j,3] - storage_energy_tracker[j,7]*1./storage_energy_tracker[j,2] #Decrease the SOC --> SOCnew = SOCold + (PowerUsed X Timestep)/MaxEnergy
                #storage_energy_tracker[j,3] = max.(storage_energy_tracker[j,3],0) #Ensure that the SOC does not drop below 0 REMOVED BECAUSE THE ABOVE CODE SHOULD ALREADY PREVENT THIS
            elseif storage_energy_tracker[j,7] < 0
                #Charge
                storage_energy_tracker[j,3] = storage_energy_tracker[j,3] - storage_energy_tracker[j,7]*1./storage_energy_tracker[j,2] #Increase the SOC (note the the value of storage_energy_tracker[j,7] should be negative)
                #storage_energy_tracker[j,3] = min.(storage_energy_tracker[j,3],1) #Ensure that the SOC does not go above 1 REMOVED BECAUSE THE ABOVE CODE SHOULD ALREADY PREVENT THIS
            else # = 0 (Hasn't been used for charging or discharging)
                #Update SOC because of losses
                #This assumes that losses only occur if no other charging/discharging occurs to a device
                storage_energy_tracker[j,3] = (1-0.01/24)*storage_energy_tracker[j,3] #This assumes about 1% loss every day
                storage_energy_tracker[j,3] = max.(storage_energy_tracker[j,3],0) #Don't let losses drop below 0
            end
        end

        #Now update the charge/discharge power columns with the new values
        storage_energy_tracker[:,5] = min.(storage_energy_tracker[:,2] .* storage_energy_tracker[:,3], storage_energy_tracker[:,1])
        storage_energy_tracker[:,6] = min.(storage_energy_tracker[:,2] .* (1 - storage_energy_tracker[:,3]), storage_energy_tracker[:,1])
        #########################################################

        #########################################################
        #Update DR here
        #Adjust payback time of any outstanding DR
        #Remove any unused rows
        DR_energy_tracker[:,2] -= 1

        temp = zeros(0,size(DR_energy_tracker,2))
        for j in 1:size(DR_energy_tracker,1)

            if (DR_energy_tracker[j,1] == 0) | (DR_energy_tracker[j,2] < 0) #If the energy/power to serve is zero or if the time has expired
                #We want to get rid of it, so don't include it in the temporary variable
            elseif DR_energy_tracker[j,1] > 0
                temp = [temp; DR_energy_tracker[j,:].']
            else
                #There shouldn't be anything "else"
            end
        end
        DR_energy_tracker = temp
        #########################################################
    println("Iteration ",i)

        #TODO: Save data better as this is inefficient to overwrite it each step
        output_data[i,1:6] = solution_vector.'
        output_data[i,7] = total_load[i]
        output_data[i,8] = sum(DR_energy_tracker[:,1])

            if size(DR_energy_tracker,1) > 0
                output_data[i,9] = minimum(DR_energy_tracker[:,2]) + 1 #Plus one is to undo the -1 performed in the DR adjustment above
            else
                output_data[i,9] = -1 + 1
            end

        output_data[i,10] = gen_upper_limits[1]
        output_data[i,11] = sum(generator_state_vector)

    end

    #Save data
    UnservedLoadData[:,MCI] = output_data[:,4]
    AvailableGenCap[:,MCI] = output_data[:,10]
    DR_Injected[:,MCI] = output_data[:,3]
    RunTime[MCI] = toc()
    #Old Code
    # df = DataFrame(:GenerationUsed=>output_data[:,1],:StorageUsed=>output_data[:,2], :DR_Injected=>output_data[:,3],
    # :UnservedLoad=>output_data[:,4],:StorageCharged=>output_data[:,5],:DR_Repayed=>output_data[:,6],
    # :Total_Load=>output_data[:,7], :DR_energy_tracker_total_energy=>output_data[:,8], :DR_energy_tracker_min_time_left=>output_data[:,9],
    # :AvailGenCap=>output_data[:,10], :FractionAvailableGenerators=>output_data[:,11]./size(gen_distributions_sequential,1))
    # CSV.write("C:/Users/aklem/Documents/GitHub/ResourceAdequacy/output.csv",df)

    println("MC Iteration", MCI)
end

###############################################################################
#Statistical Analysis
sum(UnservedLoadData,1)
UnservedHours = zeros(size(UnservedLoadData,2))
for i in 1:size(UnservedLoadData,2)
    temp = find(UnservedLoadData[:,i])
    UnservedHours[i] = length(temp)
end

###############################################################################


save("OutputData.jld","UnservedLoadData",UnservedLoadData)
save("OutputData.jld","AvailableGenCap",AvailableGenCap)
save("OutputData.jld","DR_injected",DR_Injected)

# #Test Code
# r = rand(3, 3, 3)
# aa = save("data.jld", "one", r)
# a = load("data.jld")["one"]
#
# jldopen("data.jld", "w") do file
#     write(file, "r", r)  # alternatively, say "@write file A"
# end
