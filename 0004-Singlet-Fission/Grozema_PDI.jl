using QuantumDynamics
using Plots
using LinearAlgebra


# Atomic units to SI units conversion
const thz2au = 0.0001519828500716

const invcm2au = 4.55633e-6
const mev2invcm = 8.066
const mev2au=mev2invcm*invcm2au

const eV2Ha=1/27.2114079527

const nm2au = 18.897

const au2fs = 0.02418884254

function pdi()
    println("Entering PDI routine...")
    Es = 2.22 
    Ec = 3.11 
    Et =  1.05 
    t2e =   -0.42 / 1000
    thh = -125.0  / 1000
    tll =  145.0  / 1000
    thl = -125.0  / 1000

    N = 5

    H0 = Matrix{ComplexF64}([
        Es 0 tll -thh t2e
        0 Es -thh tll t2e
        tll -thh Ec 0 √(3/2)*thl
        -thh tll 0 Ec √(3/2)*thl
        t2e t2e √(3/2)*thl √(3/2)*thl 2*Et
    ]) * eV2Ha
    println("H0 (Ha) = ", H0)

    Jw = SpectralDensities.DrudeLorentz(λ=135E-3*eV2Ha, γ=150E-3*eV2Ha, Δs=1.0)
    println("Jw = ", Jw)

    nsteps = 2000
    ρ0 = Matrix{ComplexF64}(zeros(N, N))
    ρ0[1, 1] = 1.0

    #T = 50.0:50.0:500.0

    β = 1 / (300 * 3.16683e-6) # JMF what units?

    dt = 0.05/au2fs

    d = 1 # JMF: what this?

    ωs, cs = SpectralDensities.discretize(Jw, 100)

    hb = Solvents.HarmonicBath(β, ωs, cs, [1.0, 2.0, 3.0, 4.0, 5.0].*d, 1000)
    
    #println("hb = ", hb)

    println("Propagating...")
    tc, ρc = QCPI.propagate( 
        ; Hamiltonian=H0, 
        Jw, 
        solvent=hb, 
        ρ0, 
        classical_dt=dt, dt, 
        ntimes=nsteps, 
        kmax=1, 
        svec=[1.0 2.0 3.0 4.0 5.0].*d, 
        extraargs=QuAPI.QuAPIArgs(), 
        path_integral_routine=QuAPI.propagate)
    println("Propagation done.")

    plot(tc.*au2fs, real.(ρc[:, 1, 1]), label="S1S0")
    plot!(tc.*au2fs, real.(ρc[:, 3, 3]), label="CT")
    plot!(tc.*au2fs, real.(ρc[:, 5, 5]), label="TT")

    savefig(".png")

end

pdi()

