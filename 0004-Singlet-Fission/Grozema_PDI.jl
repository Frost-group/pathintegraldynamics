using QuantumDynamics
using Plots
using LinearAlgebra



const thz2au = 0.0001519828500716
const invcm2au = 4.55633e-6
const au2fs = 0.02418884254
const mev2invcm = 8.066
const nm2au = 18.897


function pdi()
    Es = 2.22 * 1000
    Ec = 3.11 * 1000
    Et =  1.05 * 1000
    t2e = -0.42
    thh = -125.0
    tll = 145.0
    thl = -125.0

    N = 5

    H0 = Matrix{ComplexF64}([
        Es 0 tll -thh t2e
        0 Es -thh tll t2e
        tll -thh Ec 0 √(3/2)*thl
        -thh tll 0 Ec √(3/2)*thl
        t2e t2e √(3/2)*thl √(3/2)*thl 2*Et
    ]) * mev2invcm * invcm2au

    Jw = SpectralDensities.DrudeLorentz(λ=135*mev2invcm*invcm2au, γ=150*mev2invcm*invcm2au, Δs=1.0)

    nsteps = 2000
    ρ0 = Matrix{ComplexF64}(zeros(N, N))
    ρ0[1, 1] = 1.0

    T = 50.0:50.0:500.0

    β = 1 / (T[6] * 3.16683e-6) # T = 300K

    dt = 0.05/au2fs

    d = 1

    ωs, cs = SpectralDensities.discretize(Jw, 100)

    hb = Solvents.HarmonicBath(β, ωs, cs, [1.0, 2.0, 3.0, 4.0, 5.0].*d, 1000)
    tc, ρc = QCPI.propagate(; Hamiltonian=H0, Jw, solvent=hb, ρ0, classical_dt=dt, dt, ntimes=nsteps, kmax=1, svec=[1.0 2.0 3.0 4.0 5.0].*d, extraargs=QuAPI.QuAPIArgs(), path_integral_routine=QuAPI.propagate)

    plot(tc.*au2fs, real.(ρc[:, 1, 1]), label="S1S0")
    plot!(tc.*au2fs, real.(ρc[:, 3, 3]), label="CT")
    plot!(tc.*au2fs, real.(ρc[:, 5, 5]), label="TT")

    savefig(".png")

end

pdi()
