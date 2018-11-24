gens1 = reshape(RA.DispatchableGeneratorSpec.(
    [10., 10. , 10.], [0.1, 0.1, 0.1], [0.9, 0.9, 0.9]), 3, 1)

gens2 = reshape(RA.DispatchableGeneratorSpec.(
    [10., 20. , 15., 25], [0.1, 0.1, 0.1, 0.1], [0.9, 0.9, 0.9, 0.9]), 2, 2)

stors1 = Matrix{RA.StorageDeviceSpec{Float64}}(0,1)

# Single-Node System A
singlenode_a = ResourceAdequacy.SystemModel{4,1,Hour,MW,MWh}(
    gens1, stors1, DateTime(2010,1,1):Hour(1):DateTime(2010,1,1,3),
    ones(Int, 4), ones(Int, 4), [5., 6, 7, 8], [25., 28, 27, 24])

# Single-Node System B
singlenode_b = ResourceAdequacy.SystemModel{6,1,Hour,MW,MWh}(
    gens2, stors1, DateTime(2015,6,1):Hour(1):DateTime(2015,6,1,5),
    [1,1,1,2,2,2], ones(Int, 6), [7.,8,9,9,8,7], [28.,29,30,31,32,33])

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

threenode_a = ResourceAdequacy.SystemInputStateDistribution{1,Hour,MW,MWh}(
    ["A", "B", "C"],
    gen_dists, vg,
    line_labels, line_dists,
    load
)

# Three-Node System B
gen_dists = [Generic([0., 1, 2], [.1, .3, .6]),
             Generic([0., 1, 2], [.1, .3, .6]),
             Generic([0., 1, 2], [.1, .3, .6])]
gen_dists2 = [Generic([0., 1, 2], [.3, .3, .4]),
              Generic([0., 1, 2], [.1, .3, .6]),
              Generic([0., 1, 2], [.1, .2, .7])]

line_labels = [(1,2), (2,3), (1,3)]
line_dists = [Generic([0, 1.], [.2, .8]),
              Generic([0, 1.], [.2, .8]),
              Generic([0, 1.], [.2, .8])]
line_dists2 = [Generic([0, 1.], [.4, .6]),
               Generic([0, 1.], [.3, .7]),
               Generic([0, 1.], [.2, .8])]

vg = [.2 .4; 0 0; .1 .15]
load = [.5 1.5; .5 1.5; .5 1.5]

threenode_b = ResourceAdequacy.SystemInputStateDistribution{1,Hour,MW,MWh}(
    ["Region 1", "Region 2", "Region 3"], gen_dists, vg,
    line_labels, line_dists,
    load
)


# Unit-level system definitions

generators1 = [ResourceAdequacy.DispatchableGeneratorSpec(1., 0.1, 0.5),
               ResourceAdequacy.DispatchableGeneratorSpec(2., 0.1, 0.5)]
generators2 = [ResourceAdequacy.DispatchableGeneratorSpec(1., 0.3, 0.05),
               ResourceAdequacy.DispatchableGeneratorSpec(2., 0.3, 0.01)]

storages1 = [ResourceAdequacy.StorageDeviceSpec(1., 4., 0.99, 0., 1.0)]
storages2 = [ResourceAdequacy.StorageDeviceSpec(1., 4., 0.99, 0.05, 0.5)]

vg   = [0.5, 0.7, 0.3]
load = [1. , 2  , 3]

lines = [ResourceAdequacy.LineSpec(1., 0.004, 0.04),
         ResourceAdequacy.LineSpec(1., 0.017, 0.04),
         ResourceAdequacy.LineSpec(1., 0.017, 0.04)]

# Single-Node Multi-Period System
singlenode_multiperiod =
    ResourceAdequacy.SystemModel{3,1,Hour,MW,MWh}(
        hcat(generators2), hcat(storages2),
        DateTime(2020,1,1,0):Hour(1):DateTime(2020,1,1,2), ones(Int,3), ones(Int,3),
        [0.5, 0.8, 0.2], [1.,3.,1.]
    )

# Three-Node System Set

threenode_multiperiod =
    ResourceAdequacy.SystemModel{4,1,Hour,MW,MWh}(
        ["Region A", "Region B", "Region C"],
        hcat(vcat(generators1, generators2)), [1,2,4],
        hcat(vcat(storages1, storages2)), [1,3,3],
        line_labels, hcat(lines), [1,2,3],
        DateTime(2018,10,30,0):Hour(1):DateTime(2018,10,30,3),
        ones(Int, 4), ones(Int, 4), ones(Int, 4), # Time set references
        [.8 .7 .6 .7; .6 .4 .5 .7; .7 .8 .9 .8], # VG
        [1.4 1.5 1.6 1.7; 1.5 1.6 1.7 1.6; 1.3 1.4 1.5 1.6] # Load
)
