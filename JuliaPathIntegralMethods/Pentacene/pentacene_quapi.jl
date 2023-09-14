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

    Jw = SpectralDensities.DrudeLorentz(; λ=100.0*invcm2au, γ=1/(50.0 / au2fs), Δs=1.0)
    
    fbU = Propagators.calculate_bare_propagators(; Hamiltonian=H0, dt=dt, ntimes=nsteps)
    #t, ρ = QuAPI.propagate(; fbU=fbU, Jw=[Jw], β=β, ρ0=ρ0, dt=dt, ntimes=nsteps, kmax=1)
    t, ρ = TTM.propagate(; fbU=fbU, Jw=[Jw], β=β, ρ0=ρ0, dt=dt, ntimes=nsteps, rmax=1, extraargs=QuAPI.QuAPIArgs(), path_integral_routine=QuAPI.build_augmented_propagator)
    plot!(t.*au2fs, real.(ρ[:, 1,1]))
    savefig("pentacene_quapi.png")
end

pentacene_quapi()



