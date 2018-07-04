# Single-Node System A
singlenode_a = ResourceAdequacy.SystemDistribution{1,Hour,MW}(
    [Generic([2.,3,4], [.3, .4, .3])],
    zeros(1, 10),
    Tuple{Int,Int}[],
    Generic{Float64, Float64, Vector{Float64}}[],
    [1. 1 1 1 2 2 2 2 3 3]
)

# Single-Node System B
singlenode_b = ResourceAdequacy.SystemDistribution{1,Hour,MW}(
    Generic([2.,3,4], [.001, .4, .599]),
    zeros(100),
    [ones(59); fill(2., 40); 3]
)

# Three-Node System A
gen_dists = [Generic([2., 3], [.4, .6]),
             Generic([2., 3], [.4, .6]),
             Generic([2., 3], [.4, .6])]
vg = zeros(3,5)
load = Matrix{Float64}(3,5)
load[:, 1:3] = 2.
load[:, 4:5] = 2.5
line_labels = [(1,2), (2,3), (1,3)]
line_dists = [Generic([0., 1], [.1, .9]),
              Generic([0., 1], [.3, .7]),
              Generic([0., 1], [.3, .7])]

threenode_a = ResourceAdequacy.SystemDistribution{1,Hour,MW}(
    gen_dists, vg,
    line_labels, line_dists,
    load
)

# Three-Node System B
gen_dists = [Generic([0., 1, 2], [.1, .3, .6]),
             Generic([0., 1, 2], [.1, .3, .6]),
             Generic([0., 1, 2], [.1, .3, .6])]
vg = [.2 .4; 0 0; .1 .15]
line_labels = [(1,2), (2,3), (1,3)]
line_dists = [Generic([0, 1.], [.2, .8]),
              Generic([0, 1.], [.2, .8]),
              Generic([0, 1.], [.2, .8])]
load = [.5 1.5; .5 1.5; .5 1.5]

threenode_b = ResourceAdequacy.SystemDistribution{1,Hour,MW}(
    gen_dists, vg,
    line_labels, line_dists,
    load
)


# Three-Node System Set

threenode_multiperiod = ResourceAdequacy.SystemDistributionSet{1,Hour,4,Hour,MW,Float64}(
    collect(DateTime(2018,10,30,0):Dates.Hour(1):DateTime(2018,10,30,3)),
    gen_dists, [.8 .7 .6 .7; .6 .4 .5 .7; .7 .8 .9 .8],
    line_labels, line_dists,
    [1.4 1.5 1.6 1.7; 1.5 1.6 1.7 1.6; 1.3 1.4 1.5 1.6], 1, 1
)



# Three-Node System A Chronological
gen_dists = [Generic([2., 3], [.4, .6]),
             Generic([2., 3], [.4, .6]),
             Generic([2., 3], [.4, .6])]

#Generators: [MaxCap Node/Area MTTR FOR]
gen_dists = [100 1 2 0.05;
             200 1 3 0.1;
             300 2 3 0.1;
             1000 3 2 0.08;
             800 3 3 0.05]

#State transition probabilities
# 0-0 , 0-1, 1-0, 1-1
#= TODO Remove this because we will use MTTR and FOR to calculate these
gen_state_trans_probs = [0.75 0.25 0.05 0.95;
                         0.75 0.25 0.05 0.95;
                         0.75 0.25 0.05 0.95]

generator_MTTR = [2;3;2] #Mean time to repair
generator_MTBF = [20;32;15] #Mean time between failures
=#

#Storage Parameters (Max Power, Max Energy Cap, Initial SOC, Node/Area)
storage_params = [1 4 rand(1) 1;
                         2 0.5 rand(1) 1;
                         0.5 0.5 rand(1) 2
                         1 2 rand(1) 3
                         2 1 rand(1) 3]

#Demand Response Parameters
DR_params = []


vg = zeros(3,5)
load = Matrix{Float64}(3,5)
load[:, 1:3] = 2.
load[:, 4:5] = 2.5
line_labels = [(1,2), (2,3), (1,3)]
line_dists = [Generic([0., 1], [.1, .9]),
              Generic([0., 1], [.3, .7]),
              Generic([0., 1], [.3, .7])]

threenode_a = ResourceAdequacy.SystemDistribution{1,Hour,MW}(
    gen_dists, vg,
    line_labels, line_dists,
    load
)
