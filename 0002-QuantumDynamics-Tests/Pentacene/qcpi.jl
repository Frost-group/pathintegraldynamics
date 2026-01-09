using QuantumDynamics
using Plots

const thz2au = 0.0001519828500716
const invcm2au = 4.55633e-6
const au2fs = 0.02418884254


function qcpi()
    threshold = 1e-10
    nsteps = 100
    N = 4
    β = 1 / (300 * 3.16683e-6)

    Jw = SpectralDensities.DrudeLorentz(; λ=100.0*invcm2au, γ=50.0*invcm2au, Δs=1.0)
    
    dt = 100 / au2fs

    ϵb = 600 * invcm2au

    H0 = Matrix{ComplexF64}(zeros(N, N))
    ρ0 = Matrix{ComplexF64}(zeros(N, N))
    ρ0[1, 1] = 1.0

    for i in 1:N
        if i <= N-1
            H0[i, i+1] = ϵb
        end
        if i >= 2
            H0[i, i-1] = ϵb
        end
    end


    ω, c = SpectralDensities.discretize(Jw, 100)
    hb = Solvents.HarmonicBath(β, ω, c, [1.0, 2.0, 3.0, 4.0], 1000)
    tc, ρc = QCPI.propagate(; Hamiltonian=H0, Jw, solvent=hb, ρ0, classical_dt=dt / 100, dt, ntimes=nsteps, kmax=1, svec=[1.0 2.0 3.0 4.0], extraargs=QuAPI.QuAPIArgs(), path_integral_routine=QuAPI.propagate)

    plot(tc.*au2fs, real.(ρc[:, 1,1]), label="1")
    plot!(tc.*au2fs, real.(ρc[:, 2,2]), label="2")
    plot!(tc.*au2fs, real.(ρc[:, 3,3]), label="3")
    plot!(tc.*au2fs, real.(ρc[:, 4,4]), title = "QCPI", label="4")

    savefig("QCPI_Pentacene.png")
end

qcpi()
