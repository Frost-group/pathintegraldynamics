# UInt8 Truncation error reproducing

using QuantumDynamics
using DelimitedFiles

const invcm2au = 4.55633e-6
const au2fs = 0.02418884254

function QuAPITTM(; N=20, v=50*invcm2au, reorg=350*invcm2au, cutoff=40*invcm2au)
    dt = 0.25/au2fs
    nsteps = 1000
    H0 = Matrix{ComplexF64}(zeros(N, N))
    for i in 1:N
        if i <= N-1
            H0[i, i+1] = v
        end
        if i >= 2
            H0[i, i-1] = v
        end
    end

    ρ0 = Matrix{ComplexF64}(zeros(N, N))
    ρ0[1, 1] = 1.0
    β = 1 / (300 * 3.16683e-6) # T = 300K
    svec = Matrix{Float64}(zeros(1, N))
    for i in 1:N
        svec[i] = i
    end
    
    Jw = SpectralDensities.DrudeLorentz(λ=reorg, γ=cutoff, Δs=1.0)

    fbU = Propagators.calculate_bare_propagators(; Hamiltonian=H0, dt=dt, ntimes=nsteps)

    ts, ρs = TTM.propagate(; fbU=fbU,
                            Jw=[Jw],
                            β=β,
                            ρ0=ρ0,
                            dt=dt,
                            ntimes=nsteps,
                            rmax=5,
                            svec=svec,
                            extraargs=QuAPI.QuAPIArgs(),
                            path_integral_routine=QuAPI.build_augmented_propagator)
    
    
    open("populations.txt", "w") do io
        pops = [real.(ρs[:, i, i]) for i in 1:N]
        tpops = [ts pops...]
        writedlm(io, tpops, ' ')
    end
end

    

function TEMPOTTM(; N=20, v=50*invcm2au, reorg=350*invcm2au, cutoff=40*invcm2au)
    dt = 0.25/au2fs
    nsteps = 1000

    H0 = Matrix{ComplexF64}(zeros(N, N))
    for i in 1:N
        if i <= N-1
            H0[i, i+1] = v
        end
        if i >= 2
            H0[i, i-1] = v
        end
    end

    ρ0 = Matrix{ComplexF64}(zeros(N, N))
    ρ0[1, 1] = 1.0
    β = 1 / (300 * 3.16683e-6) # T = 300K
    svec = Matrix{Float64}(zeros(1, N))
    for i in 1:N
        svec[i] = i
    end
    
    Jw = SpectralDensities.DrudeLorentz(λ=reorg, γ=cutoff, Δs=1.0)

    fbU = Propagators.calculate_bare_propagators(; Hamiltonian=H0, dt=dt, ntimes=nsteps)

    ts, ρs = TTM.propagate(; fbU=fbU,
                            Jw=[Jw],
                            β=β,
                            ρ0=ρ0,
                            dt=dt,
                            ntimes=nsteps,
                            rmax=5,
                            svec=svec,
                            extraargs=TEMPO.TEMPOArgs(),
                            path_integral_routine=TEMPO.build_augmented_propagator)
    
    
    open("populations.txt", "w") do io
        pops = [real.(ρs[:, i, i]) for i in 1:N]
        tpops = [ts pops...]
        writedlm(io, tpops, ' ')
    end
end

QuAPITTM()

#TEMPOTTM()

