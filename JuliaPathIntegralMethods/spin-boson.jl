# Gao et al. JCTC 2020


using QuantumDynamics
using Plots

function spin_boson_a()
    H = Utilities.create_tls_hamiltonian(; ϵ=0.0, Δ=1.0)
    ρ0 = [1.0+0.0im 0.0; 0.0 0.0]
    
    ξ = 0.09
    ωc = 2.5
    
    Jw = SpectralDensities.ExponentialCutoff(; ξ, ωc, n=1)
    
    dt = 0.25
    ntimes = 100
    
    β = 0.1

    barefbU = Propagators.calculate_bare_propagators(; Hamiltonian=H, dt=dt, ntimes=ntimes)
    
    t, ρs = QuAPI.propagate(; fbU=barefbU, Jw=[Jw], β, ρ0, dt, ntimes, kmax=3)

    t2, ρs2 = TEMPO.propagate(; fbU=barefbU, Jw=[Jw], β, ρ0, dt, ntimes, kmax=3)

    plot!(t, real.(ρs[:,1,1] - ρs[:,2,2]), label="QuAPI")
    plot!(t2, real.(ρs2[:,1,1] - ρs2[:,2,2]), label="TEMPO")
    
    savefig("sb_a.png")

end

function spin_boson_b()
    H = Utilities.create_tls_hamiltonian(; ϵ=0.0, Δ=1.0)
    ρ0 = [1.0+0.0im 0.0; 0.0 0.0]
    
    ξ = 0.09
    ωc = 2.5
    
    Jw = SpectralDensities.ExponentialCutoff(; ξ, ωc, n=1)
    
    dt = 0.25
    ntimes = 100
    
    β = 5

    barefbU = Propagators.calculate_bare_propagators(; Hamiltonian=H, dt=dt, ntimes=ntimes)
    
    t, ρs = QuAPI.propagate(; fbU=barefbU, Jw=[Jw], β, ρ0, dt, ntimes, kmax=3)

    t2, ρs2 = TEMPO.propagate(; fbU=barefbU, Jw=[Jw], β, ρ0, dt, ntimes, kmax=3)

    plot!(t, real.(ρs[:,1,1] - ρs[:,2,2]), label="QuAPI")
    plot!(t2, real.(ρs2[:,1,1] - ρs2[:,2,2]), label="TEMPO")
    
    savefig("sb_b.png")

end


function spin_boson_c()
    H = Utilities.create_tls_hamiltonian(; ϵ=1.0, Δ=1.0)
    ρ0 = [1.0+0.0im 0.0; 0.0 0.0]
    
    ξ = 0.1
    ωc = 1.0
    
    Jw = SpectralDensities.ExponentialCutoff(; ξ, ωc, n=1)
    
    dt = 0.25
    ntimes = 100
    
    β = 0.25

    barefbU = Propagators.calculate_bare_propagators(; Hamiltonian=H, dt=dt, ntimes=ntimes)
    
    t, ρs = QuAPI.propagate(; fbU=barefbU, Jw=[Jw], β, ρ0, dt, ntimes, kmax=3)

    t2, ρs2 = TEMPO.propagate(; fbU=barefbU, Jw=[Jw], β, ρ0, dt, ntimes, kmax=3)

    plot!(t, real.(ρs[:,1,1] - ρs[:,2,2]), label="QuAPI")
    plot!(t2, real.(ρs2[:,1,1] - ρs2[:,2,2]), label="TEMPO")
    
    savefig("sb_c.png")

end

function spin_boson_d()
    H = Utilities.create_tls_hamiltonian(; ϵ=1.0, Δ=1.0)
    ρ0 = [1.0+0.0im 0.0; 0.0 0.0]
    
    ξ = 0.1
    ωc = 2.0
    
    Jw = SpectralDensities.ExponentialCutoff(; ξ, ωc, n=1)
    
    dt = 0.25
    ntimes = 100
    
    β = 5

    barefbU = Propagators.calculate_bare_propagators(; Hamiltonian=H, dt=dt, ntimes=ntimes)
    
    t, ρs = QuAPI.propagate(; fbU=barefbU, Jw=[Jw], β, ρ0, dt, ntimes, kmax=3)

    t2, ρs2 = TEMPO.propagate(; fbU=barefbU, Jw=[Jw], β, ρ0, dt, ntimes, kmax=3)

    plot!(t, real.(ρs[:,1,1] - ρs[:,2,2]), label="QuAPI")
    plot!(t2, real.(ρs2[:,1,1] - ρs2[:,2,2]), label="TEMPO")
    
    savefig("sb_d.png")

end

#spin_boson_a()
#spin_boson_b()
#spin_boson_c()
spin_boson_d()
