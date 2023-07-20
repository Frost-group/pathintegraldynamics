using QuantumDynamics
using Plots


function quapi()
    H0 = Utilities.create_tls_hamiltonian(; ϵ=0.0, Δ=2.0)
    ρ0 = [1.0+0.0im 0.0; 0.0 0.0]
    β=1.0
    dt=0.25
    ωc=1.0
    ξ=2.0

    ntimes = 600

    Jw = SpectralDensities.ExponentialCutoff(; ξ=ξ, ωc=ωc, n=1)

    fbU = Propagators.calculate_bare_propagators(; Hamiltonian=H0, dt=dt, ntimes=ntimes)

    t, ρs = QuAPI.propagate(; fbU=fbU, Jw=[Jw], β=β, ρ0=ρ0, dt=dt, ntimes=ntimes, kmax=1, svec=[0.0 2.0])

    f = open("quapi-smatpi-pathsum.txt")
    lines = readlines(f)
    
    t_smatpi = []
    ρ_smatpi = []

    for line in lines
        tokens = split(line)
        push!(t_smatpi, parse(Float64, tokens[1]))
        push!(ρ_smatpi, parse(Float64, tokens[2]))
    end

    plot!(t, real.(ρs[:,1,1]), label="QuAPI")
    plot!(t_smatpi, ρ_smatpi, label="SMatPI QuAPI")
    
    savefig("quapi-smatpi.png")
end


quapi()

