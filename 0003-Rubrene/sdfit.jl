# Testing a makeshift Spectral Density fit


using Plots

using QuantumDynamics

const thz2au = 0.0001519828500716
const invcm2au = 4.55633e-6
const au2fs = 0.02418884254
const mev2invcm = 8.066

struct fitsd <: SpectralDensities.ContinuousSpectralDensity
    ωs :: Vector{Float64}
    jws :: Vector{Float64}
end

function evaluate(sd::fitsd, ω::Real)
    ωs = sd.ωs
    jws = sd.jws
    for i in 1:(size(ωs)[1]-1)
        if ω > ωs[i] && ω < ωs[i+1]
            return jws[i] + ((jws[i+1]-jws[i])/(ωs[i+1] - ωs[i]))*(ω - ωs[i])
        end
    end
    return 0.0
end

function sdfit()
    ωp = [57.8, 59.6, 89.0, 107.3, 139.1, 639.1, 1011.2, 1344.7, 1593.3] .* invcm2au
    ωpg0p = [-1.7, 1.4, 1.6, -0.14, -2.3, -7.5, -3.6, 19.8, -42.0] .* mev2invcm * invcm2au
    g0p = ωpg0p ./ ωp
    jws = ((g0p.^(2)) ./ ωp).*(π/2)    
    
    Jw = fitsd(ωp, jws)
    

    plot!(ωp, jws)
    ω = 0:0.00001:0.007
    plot!(ω, evaluate.(Ref(Jw), ω), label="polyfit")
    savefig("specdens.png")
end

sdfit()


