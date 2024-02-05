# Simple script to execute HEOM and QuAPI-TTM for molecular trimer Singlet Fission simulations.

using QuantumDynamics
using Plots
using LinearAlgebra
using DelimitedFiles

const thz2au = 0.0001519828500716
const invcm2au = 4.55633e-6
const au2fs = 0.02418884254
const mev2invcm = 8.066
const mev2au = mev2invcm * invcm2au
const nm2au = 18.897


"""
SFtrimerHEOM

Sets up and runs HEOM calculation output- plot of populations of all sites.

Ett, Ext, Ect - site energies
Vtx, Vtc, Vxx, Vhomo, Vlumo, Vcc - Couplings
reorg - reorganization energy
cutoff - cutoff frequency
dt - timestep
nsteps - number of propagation steps
Lmax - Hierarchy length
K - number of Matsubara modes

"""

function SFtrimerHEOM(Ett, Ext, Ect, Vtx, Vtc, Vxx, Vhomo, Vlumo, Vcc, reorg, cutoff; dt=0.25/au2fs, nsteps=4000, L=3, K=2)

    N = 9

    H0 = Matrix{ComplexF64}([    ## Check this pls
        Ett 0.0 -Vtx Vtx 0.0 Vtc Vtc 0.0 0.0
        0.0 Ett 0.0 Vtx -Vtx 0.0 0.0 Vtc Vtc
        -Vtx 0.0 Ext Vxx 0.0 -Vhomo Vlumo 0.0 0.0
        Vtx Vtx Vxx Ext Vxx -Vlumo Vhomo Vhomo -Vlumo
        0.0 -Vtx 0.0 Vxx Ext 0.0 0.0 Vlumo -Vhomo
        Vtc 0.0 -Vhomo -Vlumo 0.0 Ect -Vcc -Vcc 0.0
        Vtc 0.0 Vlumo Vhomo 0.0 -Vcc Ect 0.0 -Vcc
        0.0 Vtc 0.0 Vhomo Vlumo -Vcc 0.0 Ect -Vcc
        0.0 Vtc 0.0 -Vlumo -Vhomo 0.0 -Vcc -Vcc Ect        
    ]) * mev2au

    ρ0 = Matrix{ComplexF64}(zeros(N, N))
    ρ0[4, 4] = 1.0
    β = 1 / (300 * 3.16683e-6) # T = 300K
    svec = Matrix{Float64}(zeros(1, N))
    for i in 1:N
        svec[i] = i
    end

    λs = repeat([reorg], N)
    γs = repeat([cutoff], N)
    
    JwH = Vector{SpectralDensities.DrudeLorentz}()
    sys_ops = Vector{Matrix{ComplexF64}}()
    for (j, (λ, γ)) in enumerate(zip(λs, γs))
       push!(JwH, SpectralDensities.DrudeLorentz(; λ, γ, Δs=1.0))
       op = zeros(N, N)
       op[j, j] = 1.0
       push!(sys_ops, op)
    end 
    
    times_HEOM, ρs = HEOM.propagate(;
                                    Hamiltonian=H0,
                                    ρ0,
                                    β,
                                    dt,
                                    ntimes=nsteps,
                                    Jw=JwH,
                                    sys_ops=sys_ops,
                                    num_modes=K,
                                    Lmax=L)
    
    plot(times_HEOM.*au2fs, (real.(ρs[:, 1, 1]) + real.(ρs[:, 2, 2])), title="HEOM", label="TT")
    plot!(times_HEOM.*au2fs, (real.(ρs[:, 3, 3]) + real.(ρs[:, 4, 4]) + real.(ρs[:, 5, 5])), title="HEOM", label="XT")
    plot!(times_HEOM.*au2fs, (real.(ρs[:, 6, 6]) + real.(ρs[:, 7, 7]) + real.(ρs[:, 8, 8]) + real.(ρs[:, 9, 9])), title="HEOM", label="CT")


    
    savefig("Trimer-Rubrene-HEOM-populations.png")

end


"""
SFtrimerTTM

Sets up and runs HEOM calculation output- plot of populations of all sites.

Ett, Ext, Ect - site energies
Vtx, Vtc, Vxx, Vhomo, Vlumo, Vcc - Couplings
reorg - reorganization energy
cutoff - cutoff frequency
dt - timestep
nsteps - number of propagation steps
rmax - entanglement memory length
"""

function SFtrimerTTM(Ett, Ext, Ect, Vtx, Vtc, Vxx, Vhomo, Vlumo, Vcc, reorg, cutoff; dt=0.25/au2fs, nsteps=4000, rmax=10)

    N = 9

    H0 = Matrix{ComplexF64}([    ## Check this pls
        Ett 0.0 -Vtx Vtx 0.0 Vtc Vtc 0.0 0.0
        0.0 Ett 0.0 Vtx -Vtx 0.0 0.0 Vtc Vtc
        -Vtx 0.0 Ext Vxx 0.0 -Vhomo Vlumo 0.0 0.0
        Vtx Vtx Vxx Ext Vxx -Vlumo Vhomo Vhomo -Vlumo
        0.0 -Vtx 0.0 Vxx Ext 0.0 0.0 Vlumo -Vhomo
        Vtc 0.0 -Vhomo -Vlumo 0.0 Ect -Vcc -Vcc 0.0
        Vtc 0.0 Vlumo Vhomo 0.0 -Vcc Ect 0.0 -Vcc
        0.0 Vtc 0.0 Vhomo Vlumo -Vcc 0.0 Ect -Vcc
        0.0 Vtc 0.0 -Vlumo -Vhomo 0.0 -Vcc -Vcc Ect        
    ]) * mev2au

    ρ0 = Matrix{ComplexF64}(zeros(N, N))
    ρ0[4, 4] = 1.0
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
                            rmax=rmax,
                            svec=svec,
                            extraargs=TEMPO.TEMPOArgs(),
                            path_integral_routine=TEMPO.build_augmented_propagator)
    

    plot(ts.*au2fs, (real.(ρs[:, 1, 1]) + real.(ρs[:, 2, 2])), title="TTM", label="TT")
    plot!(ts.*au2fs, (real.(ρs[:, 3, 3]) + real.(ρs[:, 4, 4]) + real.(ρs[:, 5, 5])), title="TTM", label="XT")
    plot!(ts.*au2fs, (real.(ρs[:, 6, 6]) + real.(ρs[:, 7, 7]) + real.(ρs[:, 8, 8]) + real.(ρs[:, 9, 9])), title="TTM", label="CT")


    
    savefig("Trimer-Rubrene-TTM-populations.png")

end


"""
SFdimerHEOM

Sets up and runs HEOM calculation output- plot of populations of all sites.

Ett, Ext, Ect - site energies
Vtx, Vtc, Vxx, Vhomo, Vlumo, Vcc - Couplings
reorg - reorganization energy
cutoff - cutoff frequency
dt - timestep
nsteps - number of propagation steps
Lmax - Hierarchy length
K - number of Matsubara modes

"""

function SFdimerHEOM(Ett, Ext, Ect, Vtx, Vtc, Vxx, Vhomo, Vlumo, Vcc, reorg, cutoff; dt=0.25/au2fs, nsteps=4000, L=3, K=2)

    N = 5

    H0 = Matrix{ComplexF64}([    ## Check this pls
        Ett -Vtx Vtx Vtc Vtc
        -Vtx Ext Vxx -Vhomo Vlumo
        Vtx Vxx Ext -Vlumo Vhomo
        Vtc -Vhomo -Vlumo Ect -Vcc
        Vtc Vlumo Vhomo -Vcc Ect
    ]) * mev2au

    ρ0 = Matrix{ComplexF64}(zeros(N, N))
    ρ0[2, 2] = 1.0
    β = 1 / (300 * 3.16683e-6) # T = 300K
    svec = Matrix{Float64}(zeros(1, N))
    for i in 1:N
        svec[i] = i
    end

    λs = repeat([reorg], N)
    γs = repeat([cutoff], N)
    
    JwH = Vector{SpectralDensities.DrudeLorentz}()
    sys_ops = Vector{Matrix{ComplexF64}}()
    for (j, (λ, γ)) in enumerate(zip(λs, γs))
       push!(JwH, SpectralDensities.DrudeLorentz(; λ, γ, Δs=1.0))
       op = zeros(N, N)
       op[j, j] = 1.0
       push!(sys_ops, op)
    end 
    
    times_HEOM, ρs = HEOM.propagate(;
                                    Hamiltonian=H0,
                                    ρ0,
                                    β,
                                    dt,
                                    ntimes=nsteps,
                                    Jw=JwH,
                                    sys_ops=sys_ops,
                                    num_modes=K,
                                    Lmax=L)
    
    plot(times_HEOM.*au2fs, (real.(ρs[:, 1, 1])), title="HEOM", label="TT")
    plot!(times_HEOM.*au2fs, (real.(ρs[:, 2, 2]) + real.(ρs[:, 3, 3])), title="HEOM", label="XT")
    plot!(times_HEOM.*au2fs, (real.(ρs[:, 4, 4]) + real.(ρs[:, 5, 5])), title="HEOM", label="CT")


    
    savefig("Dimer-Rubrene-HEOM-populations.png")
end


"""
SFdimerTTM

Sets up and runs HEOM calculation output- plot of populations of all sites.

Ett, Ext, Ect - site energies
Vtx, Vtc, Vxx, Vhomo, Vlumo, Vcc - Couplings
reorg - reorganization energy
cutoff - cutoff frequency
dt - timestep
nsteps - number of propagation steps
rmax - entanglement memory length
"""

function SFdimerTTM(Ett, Ext, Ect, Vtx, Vtc, Vxx, Vhomo, Vlumo, Vcc, reorg, cutoff; dt=0.25/au2fs, nsteps=4000, rmax=10)

    N = 5

    H0 = Matrix{ComplexF64}([    ## Check this pls
        Ett -Vtx Vtx Vtc Vtc
        -Vtx Ext Vxx -Vhomo Vlumo
        Vtx Vxx Ext -Vlumo Vhomo
        Vtc -Vhomo -Vlumo Ect -Vcc
        Vtc Vlumo Vhomo -Vcc Ect
    ]) * mev2au

    ρ0 = Matrix{ComplexF64}(zeros(N, N))
    ρ0[2, 2] = 1.0
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
                            rmax=rmax,
                            svec=svec,
                            extraargs=TEMPO.TEMPOArgs(),
                            path_integral_routine=TEMPO.build_augmented_propagator)
    

    plot(ts.*au2fs, (real.(ρs[:, 1, 1])), title="TTM", label="TT")
    plot!(ts.*au2fs, (real.(ρs[:, 2, 2]) + real.(ρs[:, 3, 3])), title="TTM", label="XT")
    plot!(ts.*au2fs, (real.(ρs[:, 4, 4]) + real.(ρs[:, 5, 5])), title="TTM", label="CT")


    
    savefig("Dimer-Rubrene-TTM-populations.png")

end

SFdimerTTM( 2.451*1000*mev2au, 2.311*1000*mev2au, 3.097*1000*mev2au, 0.0, 0.0, 0.079*1000*mev2au, -0.175*1000*mev2au, 0.086*1000*mev2au, 0.035*1000*mev2au, 0.71*1000*mev2au, 1/(100/au2fs)) 

SFtrimerTTM( 2.451*1000*mev2au, 2.311*1000*mev2au, 3.097*1000*mev2au, 0.0, 0.0, 0.079*1000*mev2au, -0.175*1000*mev2au, 0.086*1000*mev2au, 0.035*1000*mev2au, 0.71*1000*mev2au, 1/(100/au2fs)) 

SFdimerHEOM( 2.451*1000*mev2au, 2.311*1000*mev2au, 3.097*1000*mev2au, 0.0, 0.0, 0.079*1000*mev2au, -0.175*1000*mev2au, 0.086*1000*mev2au, 0.035*1000*mev2au, 0.71*1000*mev2au, 1/(100/au2fs)) 

SFtrimerHEOM( 2.451*1000*mev2au, 2.311*1000*mev2au, 3.097*1000*mev2au, 0.0, 0.0, 0.079*1000*mev2au, -0.175*1000*mev2au, 0.086*1000*mev2au, 0.035*1000*mev2au, 0.71*1000*mev2au, 1/(100/au2fs)) 

