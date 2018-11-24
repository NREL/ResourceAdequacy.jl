@testset "NonSequentialCopperPlate" begin

    result_1ab = assess(Backcast(), NonSequentialCopperplate(), Minimal(), singlenode_a)
    @test_broken LOLE(result_1a) ≈ LOLE{4,1,Hour}(0.06, 0.)
    @test_broken EUE(result_1a) ≈ EUE{4,1,Hour,MWh}(0.06, 0.)

    result_1ar = assess(REPRA(1,1), NonSequentialCopperplate(), Minimal(), singlenode_a)
    @test_broken LOLE(result_1a) ≈ LOLE{4,1,Hour}(0.06, 0.)
    @test_broken EUE(result_1a) ≈ EUE{4,1,Hour,MWh}(0.06, 0.)

    result_1bb = assess(Backcast(), NonSequentialCopperplate(), Minimal(), singlenode_b)
    @test_broken LOLE(result_1b) ≈ LOLE{6,1,Hour}(1e-5, 0.)
    @test_broken EUE(result_1b) ≈ EUE{6,1,Hour,MWh}(0.06, 0.)

    result_1br = assess(REPRA(1,1), NonSequentialCopperplate(), Minimal(), singlenode_b)
    @test_broken LOLE(result_1b) ≈ LOLE{6,1,Hour}(1e-5, 0.)
    @test_broken EUE(result_1b) ≈ EUE{6,1,Hour,MWh}(0.06, 0.)

    # TODO: Re-run tests with Temporal() and check period LOLPs/EUEs

    result_3mb = assess(Backcast(), NonSequentialCopperplate(),
                        Minimal(), threenode_multiperiod)
    @test_broken LOLE(result_3mb) ≈ LOLE{1,Hour}(0.1408, 0.)
    @test_broken EUE(result_3mb) ≈ EUE{1,Hour,MWh}(0.06, 0.)

    result_3mr = assess(REPRA(1,1), NonSequentialCopperplate(),
                        Minimal(), threenode_multiperiod)
    @test_broken LOLE(result_3mr) ≈ LOLE{1,Hour}(0.1408, 0.)
    @test_broken EUE(result_3mr) ≈ EUE{1,Hour,MWh}(0.06, 0.)

    # TODO: Re-run tests with SpatioTemporal() and check period and regional LOLPs/EUEs

end
