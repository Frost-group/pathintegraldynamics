# Simple script to execute Hierarchical Equations of Motion output for Holstein polaron models.

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
HolsteinPolaron1D

Sets up and runs HEOM calculation on 1D Holstein-Peierls chain. Output - plot of populations of all sites, MSD vs time, dMSD/dt vs time and mobility

N - number of sites
v - nearest-neighbours coupling
reorg - reorganization energy
cutoff - cutoff frequency
dt - timestep
nsteps - number of propagation steps
Lmax - Hierarchy length
K - number of Matsubara modes


"""

function HolsteinPolaron1D(N, v, reorg, cutoff; dt=0.25/au2fs, nsteps=400000, L=5, K=2)
    H0 = Matrix{ComplexF64}(zeros(N, N))
    for i in 1:N
        H0[i,i] = 0.0
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

    ts, ρs = HEOM.propagate(;
                                    Hamiltonian=H0,
                                    ρ0,
                                    β,
                                    dt,
                                    ntimes=nsteps,
                                    Jw=JwH,
                                    sys_ops=sys_ops,
                                    num_modes=K,
                                    Lmax=L)


    
    open("populations.txt", "w") do io
        pops = [real.(ρs[:, i, i]) for i in 1:N]
        tpops = [ts pops...]
        writedlm(io, tpops, ' ')
    end



end


"""
HolsteinPolaron1DTTM

Sets up and runs TTM calculation on 1D Holstein-Peierls chain. Output - plot of populations of all sites, MSD vs time, dMSD/dt vs time and mobility

N - number of sites
v - nearest-neighbours coupling
reorg - reorganization energy
cutoff - cutoff frequency
dt - timestep
nsteps - number of propagation steps
rmax - Memory length
E - External Field (DRAGONS) - enter it in eV/Å or something to avoid having to do charge multiplication?
a - Lattice spacing

"""

function HolsteinPolaron1DTTM(N, v, reorg, cutoff; dt=0.25/au2fs, nsteps=400000, rmax=10, E=0.0, a = 0.5*nm2au)
    H0 = Matrix{ComplexF64}(zeros(N, N))
    for i in 1:N
        H0[i,i] = 0.0 - i*a*E
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
                            rmax=rmax,
                            svec=svec,
                            extraargs=TEMPO.TEMPOArgs(),
                            path_integral_routine=TEMPO.build_augmented_propagator)
    
    
    open("populations.txt", "w") do io
        pops = [real.(ρs[:, i, i]) for i in 1:N]
        tpops = [ts pops...]
        writedlm(io, tpops, ' ')
    end

end
#for i in 3:7
#	HolsteinPolaron1D(10, 50*invcm2au, 161.5*invcm2au, 41*invcm2au; L=i, K=2)
#end

#HolsteinPolaron1D(10, 50*invcm2au, 161.5*invcm2au, 41*invcm2au; L=3, K=2)
HolsteinPolaron1DTTM(10, 50*invcm2au, 161.5*invcm2au, 41*invcm2au; rmax=10, E=0.1*mev2au)


