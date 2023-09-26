# Observing the dynamics in a single Rubrene unit cell.


using QuantumDynamics
using Plots

const thz2au = 0.0001519828500716
const invcm2au = 4.55633e-6
const au2fs = 0.02418884254
const mev2invcm = 8.066


function rubrene_cell()
    
    
    # Parameters from Ordejon (2017)
    
    ϵab = 4.1
    ϵac = 28.9
    ϵad = 0.0
    
    # Hamiltonian pieced together assuming interactions work the same at different unit cell sites.

    H0 = Matrix{ComplexF64}([   
            0.0 ϵab ϵac ϵad
            ϵab 0.0 ϵad ϵac
            ϵac ϵad 0.0 ϵab
            ϵad ϵac ϵab 0.0]) * mev2invcm * invcm2au
   
    # Dynamics params

    nsteps = 10000
    dt = 0.25 / au2fs
    ρ0 = Matrix{ComplexF64}(zeros(4, 4))
    ρ0[1, 1] = 1

    β = 1 / (300 * 3.16683e-6) # T = 300K
    Jw = SpectralDensities.ExponentialCutoff(; ξ=300000.0, ωc=57.8*invcm2au, n=0.0, Δs=1.0)


    fbU = Propagators.calculate_bare_propagators(; Hamiltonian=H0, dt=dt, ntimes=nsteps)
    
    t, ρ = TTM.propagate(; fbU=fbU, Jw=[Jw], β=β, ρ0=ρ0, dt=dt, ntimes=nsteps, rmax=1, extraargs=QuAPI.QuAPIArgs(), path_integral_routine=QuAPI.build_augmented_propagator)
    
    # Mobility Calculation based on DOI 10.1039/C9CP04770K 
    MSD = []
    for i in t
        s = 0.0
        for j in 1:4
            s += ρ[:, j, j]*((j - 5)^2)  
        end
        push!(MSD, s)
    end

    plot!(t, MSD)
    
    savefig("rubrene_unit_cell_mobility.png") 
    
end


rubrene_cell()
