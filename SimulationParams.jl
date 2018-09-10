struct SimulationParams

    solarshare::Float64
    windshare::Float64
    dravailable::Float64
    reservemargin::Float64
    nsamples::Int
    resultspath::String
    firmgencapacity::Float64

    function SimulationParams(ss::Float64, ws::Float64, dra::Float64, prm::Float64, ns::Int, rp::String, fgc::Float64)
        @assert ss >= 0.
        @assert ws >= 0.
        @assert ss + ws <= 1.
        @assert 0. <= dra <= 1.
        @assert 0. <= prm <= 1.
        @assert ns > 0
        @assert fgc >= 0.
        new(ss, ws, dra, prm, ns, rp, fgc)
    end

end


