# Observing the dynamics in a single Rubrene unit cell.


using QuantumDynamics
using Plots
using LinearAlgebra

const thz2au = 0.0001519828500716
const invcm2au = 4.55633e-6
const au2fs = 0.02418884254
const mev2invcm = 8.066

# Definitely need to replace this makeshift SD with cleaner code
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

function rubrene_1D()
    # Parameters from Ordejon (2017)
    ϵb = 134.0
    ϵ2b = -10.7
    
    # Hamiltonian with interactions with nearest 2 neighbours.   
    N = 10
    
    H0 = Matrix{ComplexF64}(zeros(N, N))

    for i in 1:N
        if i <= N-2
            H0[i, i+2] = ϵ2b
        end
        if i <= N-1
            H0[i, i+1] = ϵb
        end
        if i >= 2
            H0[i, i-1] = ϵb
        end
        if i>=3
            H0[i, i-2] = ϵ2b
        end
    end


    #print(H0)
    H0 = H0 * mev2invcm * invcm2au

    nsteps = 10000
    dt = 0.25 / au2fs
    ρ0 = Matrix{ComplexF64}(zeros(N, N))
    ρ0[5, 5] = 1

    β = 1 / (300 * 3.16683e-6) # T = 300K
    
    # Phonon modes and Couplings (Ordejon 2017)
   
    ωp = [57.8, 59.6, 89.0, 107.3, 139.1, 639.1, 1011.2, 1344.7, 1593.3] .* invcm2au
    ωpg0p = [-1.7, 1.4, 1.6, -0.14, -2.3, -7.5, -3.6, 19.8, -42.0] .* mev2invcm * invcm2au
    g0p = ωpg0p ./ ωp
    jws = ((g0p.^(2)) ./ ωp).*(π/2)
    
    Jw = SpectralDensities.ExponentialCutoff(; ξ=300000.0, ωc=57.8*invcm2au, n=0.0, Δs=1.0)
    #=
    sd = fitsd(ωp, jws)
    Jw(ω) = evaluate(sd, ω)
    ω = 0:0.00001:0.007
    plot!(ω, Jw.(ω), label="fit")
    plot!(ωp, jws)
    savefig("specdens.png")
    =#

    fbU = Propagators.calculate_bare_propagators(; Hamiltonian=H0, dt=dt, ntimes=nsteps)    
    t, ρ = TTM.propagate(; fbU=fbU, Jw=[Jw], β=β, ρ0=ρ0, dt=dt, ntimes=nsteps, rmax=1, extraargs=QuAPI.QuAPIArgs(), path_integral_routine=QuAPI.build_augmented_propagator)
    plot!(t.*au2fs, real.(ρ[:,5,5]), label="TTM")
    
    savefig("rubrene_1D.png") 
end


rubrene_1D()