using QuantumDynamics
using Plots


function quapi_ttm()
    H0 = Utilities.create_tls_hamiltonian(; ϵ=0.0, Δ=2.0)
    ρ0 = [1.0+0.0im 0.0; 0.0 0.0]
    β=1.0
    dt=0.25
    ωc=1.0
    ξ=2.0

    ntimes = 100

    Jw = SpectralDensities.ExponentialCutoff(; ξ=ξ, ωc=ωc, n=1)

    fbU = Propagators.calculate_bare_propagators(; Hamiltonian=H0, dt=dt, ntimes=ntimes)

    t, ρs = QuAPI.build_augmented_propagator_QuAPI_TTM(; fbU=fbU, Jw=[Jw], β=β, dt=dt, ntimes=ntimes)

    f = open("quapi-smatpi-pathsum.txt")
    lines = readlines(f)
    
    t_smatpi = []
    ρ_smatpi = []

    for line in lines
        tokens = split(line)
        push!(t_smatpi, parse(tokens[1], Float64))
        push!(ρ_smatpi, parse(tokens[2], Float64))
    end

    plot!(t, ρs[:,1,1], label="TTM QuAPI")
    plot!(t_smatpi, ρ_smatpi, label="SMatPI QuAPI")

end


quapi_ttm()

