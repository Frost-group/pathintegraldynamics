# Observing the dynamics in a single Rubrene unit cell.


using QuantumDynamics
using Plots

const thz2au = 0.0001519828500716
const invcm2au = 4.55633e-6
const au2fs = 0.02418884254
const mev2invcm = 8.066


function rubrene_1D()
    
    
    # Parameters from Ordejon (2017)
    
    ϵb = 134.0
    ϵ2b = -10.7
    
    # Hamiltonian with interactions with nearest 2 neighbours.   
    
    N = 10
    
    H0 = Matrix{ComplexF64}(zeros(N, N))

    for i in 1:N
        if i <= N-2
            H0[i, i+2] = ϵ2b
        end
        if i <= N-1
            H0[i, i+1] = ϵb
        end
        if i >= 2
            H0[i, i-1] = ϵb
        end
        if i>=3
            H0[i, i-2] = ϵ2b
        end
    end

    print(H0)

    H0 = H0 * mev2invcm * invcm2au

    #=H0 = Matrix{ComplexF64}([   
            0.0 ϵb ϵ2b 0.0 0.0
            ϵab 0.0 ϵb ϵ2b 0.0
            ϵac ϵad 0.0 ϵab
            ϵad ϵac ϵab 0.0]) * mev2invcm * invcm2au
    =#
    # Dynamics params

    nsteps = 10000
    dt = 0.25 / au2fs
    ρ0 = Matrix{ComplexF64}(zeros(N, N))
    ρ0[5, 5] = 1

    β = 1 / (300 * 3.16683e-6) # T = 300K

    Jw = SpectralDensities.DrudeLorentz(; λ=100.0*invcm2au, γ=50*invcm2au, Δs=1.0)   # Arbitrarily going with spectral density from Pentacene paper, will have to find params.
    
    fbU = Propagators.calculate_bare_propagators(; Hamiltonian=H0, dt=dt, ntimes=nsteps)
    
    t, ρ = TTM.propagate(; fbU=fbU, Jw=[Jw], β=β, ρ0=ρ0, dt=dt, ntimes=nsteps, rmax=1, extraargs=QuAPI.QuAPIArgs(), path_integral_routine=QuAPI.build_augmented_propagator)
    #= 
    num_points=20
    ω, c = SpectralDensities.discretize(Jw, 100)
    hb = Solvents.HarmonicBath(β, ω, c, [-2.0, -1.0, 0.0, 1.0], num_points)
 
    t_q, ρ_q = QCPI.propagate(; Hamiltonian=H0, Jw, solvent=hb, ρ0, classical_dt=dt / 100, dt, ntimes=nsteps, kmax=5, extraargs=QuAPI.QuAPIArgs(), path_integral_routine=QuAPI.propagate)
    =#
    plot!(t.*au2fs, real.(ρ[:,5,5]), label="TTM")
    

    #plot!(t_q.*au2fs, real.(ρ_q[:, 1,1]), label="QCPI")

    savefig("rubrene_1D.png") 
    
end


rubrene_1D()
