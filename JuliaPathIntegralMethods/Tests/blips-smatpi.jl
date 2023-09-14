using QuantumDynamics
using Plots


function blips()
    H0 = Utilities.create_tls_hamiltonian(; ϵ=0.0, Δ=2.0)
    ρ0 = [1.0+0.0im 0.0; 0.0 0.0]
    β=1.0
    dt=0.25
    ωc=1.0
    ξ=2.0

    ntimes = 600

    Jw = SpectralDensities.ExponentialCutoff(; ξ=ξ, ωc=ωc, n=1)

    fbU = Propagators.calculate_bare_propagators(; Hamiltonian=H0, dt=dt, ntimes=ntimes)

    t, ρs = Blip.build_augmented_propagator(; fbU=fbU, Jw=[Jw], β=β, dt=dt, ntimes=ntimes, svec=[0.0 2.0])

    f = open("blips-smatpi-pathsum.txt")
    lines = readlines(f)
    
    t_smatpi = []
    ρ_smatpi = []

    for line in lines
        tokens = split(line)
        push!(t_smatpi, parse(Float64, tokens[1]))
        push!(ρ_smatpi, parse(Float64, tokens[2]))
    end

    plot!(t, real.(ρs[:,1,1]), label="Blips")
    plot!(t_smatpi, ρ_smatpi, label="SMatPI Blips")
    
    savefig("blips-smatpi.png")
end


blips()

