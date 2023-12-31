using QuantumDynamics


function twenty_states()
    
    # Example 5 from the PathSum paper

    H = Matrix{ComplexF64}([
        0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        -1.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 -1.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 -1.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 -1.0 0.0 -1.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 -1.0 0.0 -1.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 -1.0 0.0 -1.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 -1.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 -1.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0
    ])

    dt = 0.25
    ntimes = 200

    β = 1.0
    ωc = 5.0
    ξ = 0.1

    
    ρ0 = Matrix{ComplexF64}([
        1.0+0.0im 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
    ])

    Jw = SpectralDensities.ExponentialCutoff(; ξ=ξ, ωc=ωc)

    fbU = Propagators.calculate_bare_propagators(; Hamiltonian=H, dt=dt, ntimes=ntimes)

    t, ρs = QuAPI.propagate(; fbU=fbU, Jw=[Jw], β=β, ρ0=ρ0, dt=dt, ntimes=ntimes, kmax=1)
    

end

twenty_states()
