using QuantumDynamics
using Plots


const thz2au = 0.0001519828500716
const invcm2au = 4.55633e-6
const au2fs = 0.02418884254


function pentacene_heom()

    num_modes = 2
    Lmax = 3
    threshold = 1e-10
    scaled=true
    nsteps = 1000
    dt = 1 / au2fs / nsteps

    #H0 = Utilities.create_tls_hamiltonian(; ϵ=300.0, Δ=1.0)*invcm2au
    H0 = create_nn_hamiltonian(; [300, 600].*invcm2au, couplings, periodic::Bool)
    ρ0 = [1.0+0.0im 0; 0 0]
    β = 1 / (300 * 3.16683e-6)

    Jw = SpectralDensities.DrudeLorentz(; λ=100.0*invcm2au, γ=1/(50.0 / au2fs), Δs=1.0)
    
    sys_ops = Vector{Matrix{ComplexF64}}()
    op1 = zeros(2,2)
    op1[1, 1] = 1.0

    op2 = zeros(2,2)
    op2[2, 2] = 1.0
    push!(sys_ops, op1)
    push!(sys_ops, op2)

    #t, ρ = HEOM.propagate(; Hamiltonian = H0, ρ0=ρ0, Jw=[Jw], β, ntimes=nsteps, dt, sys_ops, num_modes, Lmax, scaled, extraargs=Utilities.DiffEqArgs(; reltol=1e-6, abstol=1e-6)) 
    
    fbU = Propagators.calculate_bare_propagators(; Hamiltonian=H0, dt=dt, ntimes=nsteps)
    t, ρ = QuAPI.propagate(; fbU=fbU, Jw=[Jw], β=β, ρ0=ρ0, dt=dt, ntimes=nsteps, kmax=1)

    plot!(t, real.(ρ[:, 1,1]))
    savefig("pentacene.png")
end

pentacene_heom()



