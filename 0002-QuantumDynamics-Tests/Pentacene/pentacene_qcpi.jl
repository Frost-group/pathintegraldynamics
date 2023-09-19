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
                0.0 300.0
                300.0 0.0])*invcm2au
    ρ0 = [1.0+0.0im 0; 0 0]
    β = 1 / (300 * 3.16683e-6)

    Jw = SpectralDensities.DrudeLorentz(; λ=100.0*invcm2au, γ=50.0*invcm2au, Δs=1.0)
  
    num_points=20
    ω, c = SpectralDensities.discretize(Jw, 100)
    hb = Solvents.HarmonicBath(β, ω, c, [1.0, -1.0], num_points)
 
    t, ρ = QCPI.propagate(; Hamiltonian=H0, Jw, solvent=hb, ρ0, classical_dt=dt / 100, dt, ntimes=nsteps, kmax=5, extraargs=QuAPI.QuAPIArgs(), path_integral_routine=QuAPI.propagate)


    plot!(t.*au2fs, real.(ρ[:, 1,1]))
    savefig("pentacene_qcpi.png")
end

pentacene_quapi()



