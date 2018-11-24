@testset "TemporalResult" begin

    lolps = LOLP{1,Hour}.(rand(168)/10, rand(168)/100)
    eues = EUE{MWh,1,Hour}.(rand(168), 0.)

    # Multi-period constructor
    tstamps = DateTime(2012,4,1):Hour(1):DateTime(2012,4,7, 23)
    multiresult = ResourceAdequacy.MultiPeriodMinimalResult(
        tstamps,
        ResourceAdequacy.SinglePeriodMinimalResult{MW}.(
            lolps, eues, NonSequentialNetworkFlow(1000)),
        Backcast(),
        NonSequentialNetworkFlow(1000)
    )

    # Disallow metrics defined over different time periods
    @test_throws MethodError ResourceAdequacy.MultiPeriodMinimalResult(
        tstamps,
        ResourceAdequacy.SinglePeriodMinimalResult.(
            lolps, EUE{MWh,30,Minute}.(rand(168), 0.),
            NonSequentialNetworkFlow(1000)),
        Backcast(),
        NonSequentialNetworkFlow(1000)
    )

    # Metric constructors
    @test LOLE(multiresult) ≈ LOLE(lolps)
    @test EUE(multiresult) ≈ EUE(eues)

    @test timestamps(multiresult) == tstamps
    @test multiresult[tstamps[1]] ==
        ResourceAdequacy.SinglePeriodMinimalResult{MW}(
            lolps[1], eues[1], multiresult.simulationspec)
    @test_throws BoundsError multiresult[DateTime(2013,1,1,12)

end
