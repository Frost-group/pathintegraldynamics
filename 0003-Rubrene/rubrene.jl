using QuantumDynamics
using Plots
using LinearAlgebra


const thz2au = 0.0001519828500716
const invcm2au = 4.55633e-6
const au2fs = 0.02418884254
const mev2invcm = 8.066
const mev2au = mev2invcm * invcm2au
const nm2au = 18.897

function rubrene()
    ϵ0 = -5000*mev2invcm # This appears to not change the dynamics at all
    ϵb = 134.0
    ϵ2b = -10.7
    N = 12
    H0 = Matrix{ComplexF64}(zeros(N, N))
    for i in 1:N
        H0[i,i] = ϵ0
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
    H0 = H0 * invcm2au
    
    svec = Matrix{Float64}(zeros(1, N))
    for i in 1:N
        svec[i] = i
    end

    Jw = SpectralDensities.DrudeLorentz(λ=135*mev2au, γ=150*mev2au, Δs=1.0)

    nsteps = 10000
    ρ0 = Matrix{ComplexF64}(zeros(N, N))
    ρ0[5, 5] = 1.0
    
    T = 50.0:50.0:500.0
    
    β = 1 / (300 * 3.16683e-6) # T = 300K
    
    dt = 0.1/au2fs

    fbU = Propagators.calculate_bare_propagators(; Hamiltonian=H0, dt=dt, ntimes=nsteps)
    @time t, ρs = TTM.propagate(; fbU=fbU,
                            Jw=[Jw],
                            β=β,
                            ρ0=ρ0,
                            dt=dt,
                            ntimes=nsteps,
                            rmax=9,
                            svec=svec,
                            extraargs=TEMPO.TEMPOArgs(),
                            path_integral_routine=TEMPO.build_augmented_propagator)
    plot(t.*au2fs, real.(ρs[:, 5, 5]), label="site 5")
    plot!(t.*au2fs, real.(ρs[:, 4, 4]), label="site 4")
    plot!(t.*au2fs, real.(ρs[:, 3, 3]), label="site 3")
    plot!(t.*au2fs, real.(ρs[:, 2, 2]), label="site 2")
    plot!(t.*au2fs, real.(ρs[:, 1, 1]), label="site 1")

    savefig("rubrene-TTM.png")
        
    #=

    λs = repeat([135.0], N) * mev2au
    γs = repeat([150.0], N) * mev2au
    JwH = Vector{SpectralDensities.DrudeLorentz}()
    sys_ops = Vector{Matrix{ComplexF64}}()
    for (j, (λ, γ)) in enumerate(zip(λs, γs))
        push!(JwH, SpectralDensities.DrudeLorentz(; λ, γ, Δs=1.0))
        op = zeros(N, N)
        op[j, j] = 1.0
        push!(sys_ops, op)
    end




    times_HEOM, ρs_HEOM = HEOM.propagate(;
                                    Hamiltonian=H0,
                                    ρ0,
                                    β,
                                    dt,
                                    ntimes=nsteps,
                                    Jw=JwH,
                                    sys_ops=sys_ops,
                                    num_modes=1,
                                    Lmax=2)

    plot(times_HEOM.*au2fs, real.(ρs_HEOM[:, 5, 5]), label="site 5")
    plot!(times_HEOM.*au2fs, real.(ρs_HEOM[:, 4, 4]), label="site 4")
    plot!(times_HEOM.*au2fs, real.(ρs_HEOM[:, 3, 3]), label="site 3")
    plot!(times_HEOM.*au2fs, real.(ρs_HEOM[:, 2, 2]), label="site 2")
    plot!(times_HEOM.*au2fs, real.(ρs_HEOM[:, 1, 1]), label="site 1")
    
    savefig("rubrene-HEOM.png")
    =#
end

rubrene()



