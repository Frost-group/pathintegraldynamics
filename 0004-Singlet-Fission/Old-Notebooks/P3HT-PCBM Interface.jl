using QuantumDynamics
using Plots
using LinearAlgebra

const thz2au = 0.0001519828500716
const invcm2au = 4.55633e-6
const au2fs = 0.02418884254
const mev2invcm = 8.066
const mev2au = mev2invcm * invcm2au
const nm2au = 18.897

    H0 = Matrix{ComplexF64}([
        0.100 0.100 0.357 0.139 0.200 0.007
        0.100 0.100 0.139 0.357 0.014 0.013
        0.357 0.139 0.280 0.001 0.005 0.002
        0.139 0.357 0.001 0.230 0.019 0.165
        0.200 0.014 0.005 0.019 0.000 0.102
        0.007 0.013 0.002 0.165 0.102 0.140
    ]) * 1000 * mev2au    # Simplified Hamiltonian for 2 donor sites (the actual Hamiltonian is a larger extension and includes other couplings)

N=6

λs = repeat([150.0], 5) * mev2invcm * invcm2au
γs = repeat([180.0], 5) * mev2invcm * invcm2au
JwH = Vector{SpectralDensities.DrudeLorentz}()
sys_ops = Vector{Matrix{ComplexF64}}()
for (j, (λ, γ)) in enumerate(zip(λs, γs))
    push!(JwH, SpectralDensities.DrudeLorentz(; λ, γ, Δs=1.0))
    op = zeros(5, 5)
    op[j, j] = 1.0
    push!(sys_ops, op)
end
 
nsteps = 30000
ρ0 = Matrix{ComplexF64}(zeros(N, N))
ρ0[1, 1] = 1.0

T = 50.0:50.0:500.0

β = 1 / (300 * 3.16683e-6) # T = 300K

dt = 0.1/au2fs

@time times_HEOM, ρs_HEOM = HEOM.propagate(; Hamiltonian=H0, ρ0, β, dt, ntimes=nsteps, Jw=JwH, sys_ops=sys_ops, num_modes=2, Lmax=3)
