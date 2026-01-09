using QuantumDynamics
using Plots
using LinearAlgebra

const thz2au = 0.0001519828500716
const invcm2au = 4.55633e-6
const au2fs = 0.02418884254
const mev2invcm = 8.066
const mev2au = mev2invcm * invcm2au
const nm2au = 18.897

N = 5

    H0 = Matrix{ComplexF64}([
        4.144 0.026 0.109 -0.092 0.000
        0.026 4.070 0.090 -0.112 0.000
        0.109 0.090 4.509 0.000 0.070
        -0.092 -0.112 0.000 4.992 0.087
        0.000 0.000 0.070 0.087 4.065 
    ]) * 1000 * mev2au

show(stdout, "text/plain", real.(H0))

dt = 0.1/au2fs 

nsteps = 30000
ρ0 = Matrix{ComplexF64}(zeros(N, N))
ρ0[1, 1] = 1.0

T = 50.0:50.0:500.0

β = 1 / (T[6] * 3.16683e-6) # T = 300K

Jw2 = SpectralDensities.DrudeLorentz(λ=150*invcm2au, γ=180*invcm2au, Δs=1.0)
ω = 0:0.0001:0.07
ωs2, cs2 = SpectralDensities.discretize(Jw2, 100)
plot(ω, Jw2.(ω), xlabel="ω", ylabel="J(ω)", title="Spectral Density", legend=false)

using DelimitedFiles
open("Jw_1.txt", "w") do io
    writedlm(io, [ω  Jw2.(ω)], ' ')
end

fbU = Propagators.calculate_bare_propagators(; Hamiltonian=H0, dt, ntimes=nsteps)
tquapi, ρquapi = QuAPI.propagate(; 
                                fbU=fbU, 
                                Jw=[Jw2],
                                β=β,
                                ρ0=ρ0,
                                dt,
                                ntimes=nsteps,
                                svec=[1.0 2.0 3.0 4.0 5.0],
                                kmax=3,
                                verbose=true)

plot(tquapi.*au2fs, real.(ρquapi[:, 1, 1]), label="LE1")
plot!(tquapi.*au2fs, real.(ρquapi[:, 2, 2]), label="LE2")
#plot!(tquapi.*au2fs, real.(ρquapi[:, 3, 3]), label="CT1")
#plot!(tquapi.*au2fs, real.(ρquapi[:, 4, 4]), label="CT2")
plot!(tquapi.*au2fs, real.(ρquapi[:, 5, 5]), label="TT")

fbU = Propagators.calculate_bare_propagators(; Hamiltonian=H0, dt=dt, ntimes=nsteps)
@time t_TTM, ρ_TTM = TTM.propagate(; fbU=fbU,
                            Jw=[Jw2],
                            β=β,
                            ρ0=ρ0,
                            dt=dt,
                            ntimes=nsteps,
                            rmax=20,
                            svec=[1.0 2.0 3.0 4.0 5.0],
                            extraargs=TEMPO.TEMPOArgs(; cutoff=1e-15, maxdim=1000),
                            path_integral_routine=TEMPO.build_augmented_propagator)

plot(t_TTM.*au2fs, real.(ρ_TTM[:, 1, 1]), label="LE1")
plot!(t_TTM.*au2fs, real.(ρ_TTM[:, 2, 2]), label="LE2")
#plot!(t_TTM.*au2fs, real.(ρ_TTM[:, 3, 3]), label="CT1")
#plot!(t_TTM.*au2fs, real.(ρ_TTM[:, 4, 4]), label="CT2")
plot!(t_TTM.*au2fs, real.(ρ_TTM[:, 5, 5]), label="TT")

λs = repeat([150.0], 5) * mev2invcm * invcm2au
γs = repeat([180.0], 5) * mev2invcm * invcm2au
JwH = Vector{SpectralDensities.DrudeLorentz}()
sys_ops = Vector{Matrix{ComplexF64}}()
for (j, (λ, γ)) in enumerate(zip(λs, γs))
    push!(JwH, SpectralDensities.DrudeLorentz(; λ, γ, Δs=1.0))
    op = zeros(5, 5)
    op[j, j] = 1.0
    push!(sys_ops, op)
end
 
nsteps = 30000
ρ0 = Matrix{ComplexF64}(zeros(N, N))
ρ0[1, 1] = 1.0

T = 50.0:50.0:500.0

β = 1 / (300 * 3.16683e-6) # T = 300K

dt = 0.1/au2fs

@time times_HEOM, ρs_HEOM = HEOM.propagate(; Hamiltonian=H0, ρ0, β, dt, ntimes=nsteps, Jw=JwH, sys_ops=sys_ops, num_modes=2, Lmax=3)

plot(times_HEOM.*au2fs, real.(ρs_HEOM[:, 1, 1]), label="LE1")
plot!(times_HEOM.*au2fs, real.(ρs_HEOM[:, 2, 2]), label="LE2")
#plot!(times_HEOM.*au2fs, real.(ρs_HEOM[:, 3, 3]), label="CT1")
#plot!(times_HEOM.*au2fs, real.(ρs_HEOM[:, 4, 4]), label="CT2")
plot!(times_HEOM.*au2fs, real.(ρs_HEOM[:, 5, 5]), label="TT")

@time times_BRME, ρs_BRME = BlochRedfield.propagate(;
                                            Hamiltonian=H0,
                                            Jw=JwH,
                                            β,
                                            ρ0,
                                            dt=dt,
                                            ntimes=nsteps,
                                            sys_ops)

plot(times_BRME.*au2fs, real.(ρs_BRME[:, 1, 1]), label="LE1")
plot!(times_BRME.*au2fs, real.(ρs_BRME[:, 2, 2]), label="LE2")
#plot!(times_BRME.*au2fs, real.(ρs_BRME[:, 3, 3]), label="CT1")
#plot!(times_BRME.*au2fs, real.(ρs_BRME[:, 4, 4]), label="CT2")
plot!(times_BRME.*au2fs, real.(ρs_BRME[:, 5, 5]), label="TT")

plot(times_BRME.*au2fs, real.(ρs_BRME[:, 5, 5]), label="TT")


