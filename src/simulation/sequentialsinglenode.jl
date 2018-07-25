using MathProgBase
using Clp
using RowEchelon

#Inputs:
#Not running the assess function, just running the script
#Sequential one-node system (ignoring transmission)
timesteps = 10

region_labels = ["A","B","C"]

gen_dists = [Generic([2., 3], [.4, .6]),
             Generic([2., 3], [.4, .6]),
             Generic([2., 3], [.4, .6])]

#Generators: [MaxCap Node/Area MTTR FOR]
gen_distributions_sequential = [100 1 2 0.05;
             200 1 3 0.1;
             300 1 3 0.1;
             1000 1 2 0.08;
             800 1 3 0.05]


#Storage Parameters (Max Power, Max Energy Cap, Initial SOC, Node/Area)
storage_params = [1 4 rand(1) 1;
                         2 0.5 rand(1) 1;
                         0.5 0.5 rand(1) 1
                         1 2 rand(1) 1
                         2 1 rand(1) 1]

#Demand Response Parameters (Power, Shiftable Periods, Time periods until payback is required, Node/Area)
DR_params = [1 2 8 1
            2 4 Inf 1]


vg = zeros(1,5)
load = Matrix{Float64}(1,5)
load[:, 1:3] = 2.
load[:, 4:5] = 2.5
# line_labels = [(1,2), (2,3), (1,3)]
# line_dists = [Generic([0., 1], [.1, .9]),
#               Generic([0., 1], [.3, .7]),
#               Generic([0., 1], [.3, .7])]



#TODO: Currently hardcoded in as one T unit (e.g 1 hour, 1 minute). Fix this.

    n_gens = size(gen_distributions_sequential,1) #Number of system dispatchable generators, not necessarily equal to the number of nodes/areas
    n_storage = size(storage_params,1) #Number of energy storage devices

    lol_count = 0
    lol_sum = 0.
    lol_count_with_storage = 0


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
    system.storage_params min.(
        system.storage_params[:,2] .* system.storage_params[:,3], system.storage_params[:,1]
    ) min.(
        system.storage_params[:,2] .* (1 - system.storage_params[:,3]), system.storage_params[:,1]
    )] #Last column is min(energy X SOC, MaxPowerCap)
storage_energy_tracker = [storage_energy_tracker zeros(n_storage,1)] #The last column is to keep track of updates that will be implemented at the end of the loop. Originally, a vector was used, but with the sorting process that occurs, this can get difficult to manage
#storage_energy_tracker = [MaxPower MaxEnergy SOC Node AvailDischargeCap AvailChargeCap PowerUsedThisTimestep]
#TODO: Make storage_energy_tracker into a structure? Probably easier than a matrix
###########################################################################

###########################################################################
#Initialize DR here


#TODO: Put this into a structure?
gen_matrix = eye(n)
storage_discharge_matrix = eye(n)
DR_inject_matrix = eye(n)
unserved_load_matrix = eye(n)
storage_charge_matrix = eye(n)
DR_payback_matrix = eye(n)

#Generate the A matrix, rather the matrix of coefficients for linear optimization
#Group the elements together in the matrix i.e. [g1;g2,...;s1;s2;...;d1;d2;...]
A = [gen_matrix storage_discharge_matrix DR_inject_matrix unserved_load_matrix storage_charge_matrix DR_payback_matrix]

#Create the objective vector. Initialize the length as the width of A (the number of objective variables)
gen_cost = 0
storage_discharge_cost = 1
DR_inject_cost = 2
unserved_load_cost = 1000
storage_charge_cost = -0.5
DR_payback_cost = -1.5

c = [gen_cost*ones(n,1); storage_discharge_cost*ones(n,1);
DR_inject_cost*ones(n,1); unserved_load_cost*ones(n,1); storage_charge_cost*ones(n,1);
DR_payback_cost*ones(n,1);
]

gen_lower_limits = zeros(n,1)
storage_discharge_lower_limits = zeros(n,1)
DR_inject_lower_limits = zeros(n,1)
unserved_load_lower_limits = zeros(n,1)
storage_charge_lower_limits = zeros(n,1)
DR_payback_lower_limits = zeros(n,1)

for i in 1:timesteps
    #For our case, lower bound = upper bound, and Ax = Load Demanded, not zero
    #Thus, Ax = lb = ub = Load
    lb = LOADVECTOR #NEED TO DEFINE THIS!!!!!
    ub = lb

    #Sum all of the generation in each node, changes each time iteration
    gen_upper_limits = zeros(n,1)
    for j in 1:n
        temp_index_vector = find(gen_distributions_sequential[:,2].==j) #temporary vector of the indices of the generators at node j
        gen_upper_limits[j] = sum(gen_distributions_sequential[temp_index_vector].*generator_state_vector[temp_index_vector]) #Multiply the generator capacities by the state vector, which will either multiply by zero or one, then sum. Even though system.gen_dists is a matrix, the indices will extract the values in the first column, as Julia is column-major
    end

    #Sum all of the storage in each node, changes each time iteration
    storage_discharge_upper_limits = zeros(n,1)
    for j in 1:n
        temp_index_vector = find(storage_energy_tracker[:,4].==j) #temporary vector of the indices of the storage devices at node j
        temp_var = storage_energy_tracker[:,5] #Temporarily extract just the fifth column containing the storage_energy_tracker matrix, which contains just the max available power for this timestep for each storage device, as the max avail power may change each iteration
        storage_discharge_upper_limits[j] = sum(temp_var[temp_index_vector]) #temp_index_vector indexes into the "fifth column" of storage_energy_tracker which contains the power capacity available at this timestep, calculated from SOC and MaxEnergy
    end

    #TODO:
    DR_inject_upper_limits =


    unserved_load_upper_limtis = LOADVECTOR #Can't be any higher than the laod

    storage_charge_upper_limits = zeros(n,1)
    for j in 1:n
        temp_index_vector = find(storage_energy_tracker[:,4].==j) #temporary vector of the indices of the storage devices at node j
        temp_var = storage_energy_tracker[:,6] #Temporarily extract just the sixth column containing the storage_energy_tracker matrix, which contains just the max available power for this timestep for each storage device, as the max avail power may change each iteration
        storage_charge_upper_limits[j] = sum(temp_var[temp_index_vector]) #temp_index_vector indexes into the "sixth column" of storage_energy_tracker which contains the power capacity available at this timestep, calculated from SOC and MaxEnergy
    end

    #TODO:
    DR_payback_upper_limits =



    gen_upper_limits
    storage_discharge_upper_limits
    DR_inject_upper_limits
    unserved_load_upper_limtis
    storage_charge_upper_limits
    DR_payback_upper_limits
    #transmission_upper_limits, determined above

    l = [gen_lower_limits; storage_discharge_lower_limits; DR_inject_lower_limits;
        unserved_load_lower_limits; storage_charge_lower_limits;
        DR_payback_lower_limits]

    u = [gen_upper_limits; storage_discharge_upper_limits; DR_inject_upper_limits;
        unserved_load_upper_limits; storage_charge_upper_limits;
        DR_payback_upper_limits]


    solver = ClpSolver()
    solution = linprog(c,A,lb,ub,l,u,solver)

    #What kind of data do we want to save?

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
#########################################################

#########################################################
#Update DR here
