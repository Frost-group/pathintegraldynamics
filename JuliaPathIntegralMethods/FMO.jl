using QuantumDynamics
using Plots

function FMO_HEOM(num_modes, Lmax, scaled=true)
    
    thz2au = 0.0001519828500716
    invcm2au = 4.55633e-6
    au2fs = 0.02418884254
    
    β = 1 / (77 * 3.16683e-6)

    H = Matrix{ComplexF64}([
        12410 -87.7 5.5 -5.9 6.7 -13.7 -9.9
        -87.7 12530 30.8 8.2 0.7 11.8 4.3
        5.5 30.8 12210 -53.5 -2.2 -9.6 6.0
        -5.9 8.2 -53.5 12320 -70.7 -17.0 -63.3
        6.7 0.7 -2.2 -70.7 12480 81.1 -1.3
        -13.7 11.8 -9.6 -17.0 81.1 12630 39.7
        -9.9 4.3 6.0 -63.3 -1.3 39.7 12440
    ]) * invcm2au

    nsteps = 500
    dt = 1000 / au2fs / nsteps
    ρ0 = Matrix{ComplexF64}(zeros(7, 7))
    ρ0[1, 1] = 1

    λs = repeat([35.0], 7) * invcm2au
    γs = 1 ./ (repeat([50.0], 7) ./ au2fs)
    Jw = Vector{SpectralDensities.DrudeLorentz}()
    sys_ops = Vector{Matrix{ComplexF64}}()
    for (j, (λ, γ)) in enumerate(zip(λs, γs))
        push!(Jw, SpectralDensities.DrudeLorentz(; λ, γ, Δs=1.0))
        op = zeros(7, 7)
        op[j, j] = 1.0
        push!(sys_ops, op)
    end
    
    threshold=1e-5

    @time t, ρ = HEOM.propagate(; Hamiltonian=H, ρ0=ρ0, Jw, β, ntimes=nsteps, dt, sys_ops, num_modes, Lmax, extraargs=Utilities.DiffEqArgs(; reltol=1e-6, abstol=1e-6))
    t .*= au2fs
    #fig, ax = new_figure("full")
    for j = 1:7
        plot!(t, real.(ρ[:, j, j]), label="Site $(j)", xlabel="t", ylabel="P(t)")
    end
    savefig("heom_fmo.png")
end


function FMO_QuAPI() 
    thz2au = 0.0001519828500716
    invcm2au = 4.55633e-6
    au2fs = 0.02418884254
    
    β = 1 / (77 * 3.16683e-6)

    H = Matrix{ComplexF64}([
        12410 -87.7 5.5 -5.9 6.7 -13.7 -9.9
        -87.7 12530 30.8 8.2 0.7 11.8 4.3
        5.5 30.8 12210 -53.5 -2.2 -9.6 6.0
        -5.9 8.2 -53.5 12320 -70.7 -17.0 -63.3
        6.7 0.7 -2.2 -70.7 12480 81.1 -1.3
        -13.7 11.8 -9.6 -17.0 81.1 12630 39.7
        -9.9 4.3 6.0 -63.3 -1.3 39.7 12440
    ]) * invcm2au

    nsteps = 500
    dt = 1000 / au2fs / nsteps
    ρ0 = Matrix{ComplexF64}(zeros(7, 7))
    ρ0[1, 1] = 1

    λs = repeat([35.0], 7) * invcm2au
    γs = 1 ./ (repeat([50.0], 7) ./ au2fs)
    Jw = Vector{SpectralDensities.DrudeLorentz}()
    #sys_ops = Vector{Matrix{ComplexF64}}()
    for (j, (λ, γ)) in enumerate(zip(λs, γs))
        push!(Jw, SpectralDensities.DrudeLorentz(; λ, γ, Δs=1.0))
        #op = zeros(7, 7)
        #op[j, j] = 1.0
        #push!(sys_ops, op)
    end

    threshold=1e-5

    #Jw = SpectralDensities.ExponentialCutoff(; ξ=0.16, ωc=7.5)
    
    
    barefbU = Propagators.calculate_bare_propagators(; Hamiltonian=H, dt=dt, ntimes=nsteps)

    svec = [1.0 -1.0 ; 1.0 -1.0; 1.0 -1.0; 1.0 -1.0 ; 1.0 -1.0 ; 1.0 -1.0; 1.0 -1.0]
    
    @time times, ρ = TTM.propagate(; fbU=barefbU, ρ0=ρ0, Jw=Jw, β, ntimes=nsteps, dt, svec, rmax=1, extraargs=QuAPI.QuAPIArgs(), path_integral_routine=QuAPI.build_augmented_propagator)

    for j = 1:7
        plot!(times, real.(ρ[:, j, j]), label="Site $(j)", xlabel="t", ylabel="P(t)")
    end
    savefig("fmo_quapi_ttm.png")

end


#FMO_HEOM(2,3)

FMO_QuAPI()

