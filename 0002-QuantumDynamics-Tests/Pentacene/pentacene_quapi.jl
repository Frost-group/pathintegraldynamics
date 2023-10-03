using QuantumDynamics
using Plots


const thz2au = 0.0001519828500716
const invcm2au = 4.55633e-6
const au2fs = 0.02418884254


function pentacene_quapi()
    threshold = 1e-10
    nsteps = 1000
    dt = 1000 / au2fs / nsteps

    H0 = Matrix{ComplexF64}([
                0.0 600.0
                600.0 0.0])*invcm2au
    ρ0 = [1.0+0.0im 0; 0 0]
    β = 1 / (300 * 3.16683e-6)

    Jw = SpectralDensities.DrudeLorentz(; λ=100.0*invcm2au, γ=50.0*invcm2au, Δs=1.0)
    
    

    fbU = Propagators.calculate_bare_propagators(; Hamiltonian=H0, dt=dt, ntimes=nsteps)
    
    t, ρ = TTM.propagate(; fbU=fbU, Jw=[Jw], β=β, ρ0=ρ0, dt=dt, ntimes=nsteps, rmax=1, extraargs=QuAPI.QuAPIArgs(), path_integral_routine=QuAPI.build_augmented_propagator)
    
    for j in 1:N
        #push!(JwD, SpectralDensities.DrudeLorentz(; λ, γ, Δs=1.0))
        op = zeros(N, N)
        op[j, j] = 1.0
        push!(sys_ops, op)
    end
    

    @time th, ρh = HEOM.propagate(; Hamiltonian=H0, ρ0=ρ0, Jw=[Jw], β, ntimes=nsteps, dt, sys_ops, num_modes, Lmax, scaled, threshold, extraargs=Utilities.DiffEqArgs(; reltol=1e-6, abstol=1e-6))
    th .*= au2fs
    

    plot!(t.*au2fs, real.(ρ[:, 1,1]), label="TTM")
    plot!(th, real.(ρ[:, 1,1]), label="HEOM")
    savefig("pentacene_quapi_HEOM.png")
end

pentacene_quapi()



