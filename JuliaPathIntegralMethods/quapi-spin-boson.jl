using QuantumDynamics
using Plots

function quapi_spin_boson()
    H = Utilities.create_tls_hamiltonian(; ϵ=0.0, Δ=2.0)
    ρ0 = [1.0+0.0im 0.0; 0.0 0.0]
    
    ξ = 0.09
    ωc = 2.5
    
    Jw = SpectralDensities.ExponentialCutoff(; ξ, ωc, n=1)
    
    dt = 0.25
    ntimes = 100
    
    β = 5

    barefbU = Propagators.calculate_bare_propagators(; Hamiltonian=H, dt=dt, ntimes=ntimes)
    
    t, ρs = QuAPI.propagate(; fbU=barefbU, Jw=[Jw], β, ρ0, dt, ntimes, kmax=1)
    
    plot(t, real.(ρs[:,1,1] - ρs[:,2,2]))
    
    savefig("quapi-sb.png")

end

quapi_spin_boson()
