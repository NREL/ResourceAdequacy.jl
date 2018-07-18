@testset "SequentialNetworkFlow" begin

    @test true

    #Sequential system is defined in systems.jl
    #Do we need to include ResultSpec somehow?
    result_sequential_network = assess(SequentialNetworkFlow(100_000),threenode_seq_a)

end
