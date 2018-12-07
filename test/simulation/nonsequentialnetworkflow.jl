@testset "NonSequentialNetworkFlow" begin

    # TODO: Decide on how to programmatically test MC estimates vs ground truth
    @test_broken false

    println("Single-region system A")
    println("Theoretical:")
    println("LOLE = 0.355")
    println("EUE = 1.59")

    println("Minimal:")
    result_1ab =
        assess(Backcast(), NonSequentialNetworkFlow(100_000), Minimal(), singlenode_a)
    println(LOLE(result_1ab))
    println(EUE(result_1ab))

    println("Spatial:")
    result_1ab =
        assess(Backcast(), NonSequentialNetworkFlow(100_000), Spatial(), singlenode_a)
    println(LOLE(result_1ab), " ", LOLE(result_1ab, "Region"))
    println(EUE(result_1ab), " ", EUE(result_1ab, "Region"))

    # TODO: Check per-period results as well
    println("Temporal:")
    result_1ab =
        assess(Backcast(), NonSequentialNetworkFlow(100_000), Temporal(), singlenode_a)
    println(LOLE(result_1ab))
    println(EUE(result_1ab))

    println("\nSingle-region system B")
    println("Theoretical:")
    println("LOLE = 0.96")
    println("EUE = 7.11")

    println("Minimal:")
    result_1bb =
        assess(Backcast(), NonSequentialNetworkFlow(1_000_000), Minimal(), singlenode_b)
    println(LOLE(result_1bb))
    println(EUE(result_1bb))

    println("Spatial:")
    result_1bb =
        assess(Backcast(), NonSequentialNetworkFlow(100_000), Spatial(), singlenode_b)
    println(LOLE(result_1bb), " ", LOLE(result_1bb, "Region"))
    println(EUE(result_1bb), " ", EUE(result_1bb, "Region"))

    # TODO: Check per-period results as well
    println("Temporal:")
    result_1bb =
        assess(Backcast(), NonSequentialNetworkFlow(100_000), Temporal(), singlenode_b)
    println(LOLE(result_1bb))
    println(EUE(result_1bb))

    println("Three-region system")
    #TODO: Network case is tractable, calculate true values
    println("Theoretical:")
    println("LOLE = _")
    println("EUE = _")

    println("Minimal, Backcast:")
    result_3mb_minimal = assess(Backcast(), NonSequentialNetworkFlow(100_000),
                                Minimal(), threenode)
    println(LOLE(result_3mb_minimal))
    println(EUE(result_3mb_minimal))

    println("Minimal, REPRA(1,1):")
    result_3mr_minimal = assess(REPRA(1,1), NonSequentialNetworkFlow(100_000),
                                Minimal(), threenode)
    println(LOLE(result_3mr_minimal))
    println(EUE(result_3mr_minimal))

    # println("Network, Backcast:")
    # result_3mb_network = assess(Backcast(), NonSequentialNetworkFlow(100_000),
    #                             NetworkResult(failuresonly=true),
    #                             threenode)
    # println(LOLE(result_3mb_network))
    # println(EUE(result_3mb_network))

    # println("Network, REPRA(1,1): ")
    # result_3mr_network = assess(REPRA(1,1), NonSequentialNetworkFlow(100_000),
    #                             NetworkResult(failuresonly=true),
    #                             threenode)
    # println(LOLE(result_3mr_network))
    # println(EUE(result_3mr_network))
    # println()

end
