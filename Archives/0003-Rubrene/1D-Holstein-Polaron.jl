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

function HolsteinPolaron1D(N, v, reorg, cutoff; dt=0.25/au2fs, nsteps=10000, L=5, K=3)
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

    nsteps = 40000
    ρ0 = Matrix{ComplexF64}(zeros(N, N))
    ρ0[1, 1] = 1.0
    β = 1 / (300 * 3.16683e-6) # T = 300K
    dt = 0.25/au2fs
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

    for i in 1:N
        plot!(times_HEOM.*au2fs, real.(ρs[:, i, i]), label="site $i", xscale=:log10)
    end

    savefig("rubrene-HEOM-populations.png")

    MSD = []
    for i in 1:nsteps
        s = 0.0
        for j in 1:N
            s += real(ρs[i, j, j])*(j)^2
        end
        push!(MSD, s)
    end
    
    plot((times_HEOM[1:nsteps]).*au2fs, MSD, xscale=:log10)
    
    savefig("rubrene-HEOM-MSD.png")

    dMSD_dt = [(MSD[i] - MSD[i-1])/(dt*au2fs) for i in 2:40000]
    plot((times_HEOM[2:40000]).*au2fs, dMSD_dt, xscale=:log10)
    
    savefig("rubrene-HEOM-dMSD_dt.png")
end
