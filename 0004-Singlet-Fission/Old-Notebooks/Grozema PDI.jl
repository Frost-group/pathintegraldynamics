using QuantumDynamics
using Plots
using LinearAlgebra

const thz2au = 0.0001519828500716
const invcm2au = 4.55633e-6
const au2fs = 0.02418884254
const mev2invcm = 8.066
const mev2au = mev2invcm * invcm2au
const nm2au = 18.897

Es = 2.22 * 1000
Ec = 3.11 * 1000
Et =  1.05 * 1000
t2e = -0.42
thh = -125.0
tll = 145.0
thl = -125.0

N = 5

    H0 = Matrix{ComplexF64}([
        Es 0 tll -thh t2e
        0 Es -thh tll t2e
        tll -thh Ec 0 √(3/2)*thl
        -thh tll 0 Ec √(3/2)*thl
        t2e t2e √(3/2)*thl √(3/2)*thl 2*Et
    ]) * mev2au

show(stdout, "text/plain", real.(H0))

Jw = SpectralDensities.DrudeLorentz(λ=135*mev2au, γ=150*mev2au, Δs=1.0)

ω = 0:0.0001:0.5
plot(ω, Jw.(ω), xlabel="ω", ylabel="J(ω)", title="Spectral Density", legend=false)

using DelimitedFiles
writedlm("Jw_1.txt", zip(ω, Jw.(ω)))

nsteps = 200
ρ0 = Matrix{ComplexF64}(zeros(N, N))
ρ0[1, 1] = 1.0

T = 50.0:50.0:500.0

β = 1 / (300 * 3.16683e-6) # T = 300K

dt = 0.5/au2fs 

#ωs, cs = SpectralDensities.discretize(Jw, 100)

# ## QCPI Run

# #hb = Solvents.HarmonicBath(β, ωp, ωpg0p, [1.0, 2.0, 3.0, 4.0, 5.0].*d, 1000)
# hb = Solvents.HarmonicBath(β, ωs, cs, [1.0, 2.0, 3.0, 4.0, 5.0], 1000)
# tc, ρc = QCPI.propagate(; Hamiltonian=H0, Jw=Jw, solvent=hb, ρ0, classical_dt=dt / 100, dt, ntimes=nsteps, kmax=3, svec=[1.0 2.0 3.0 4.0 5.0], extraargs=QuAPI.QuAPIArgs(), path_integral_routine=QuAPI.propagate, verbose=true)

# plot(tc.*au2fs, real.(ρc[:, 1, 1]), label="S1S0")
# plot!(tc.*au2fs, real.(ρc[:, 3, 3]), label="CT")
# plot!(tc.*au2fs, real.(ρc[:, 5, 5]), label="TT")

# fbU = Propagators.calculate_bare_propagators(; Hamiltonian=H0, dt, ntimes=nsteps)
# tquapi, ρquapi = QuAPI.propagate(; 
#                                 fbU=fbU, 
#                                 Jw=[Jw],
#                                 β=β,
#                                 ρ0=ρ0,
#                                 dt,
#                                 ntimes=nsteps,
#                                 svec=[1.0 2.0 3.0 4.0 5.0],
#                                 kmax=3,
#                                 verbose=true)

#fbU = Propagators.calculate_bare_propagators(; Hamiltonian=H0, dt=dt, ntimes=nsteps)
#=@time t_TTM, ρ_TTM = TTM.propagate(; fbU=fbU,
                            Jw=[Jw],
                            β=β,
                            ρ0=ρ0,
                            dt=dt,
                            ntimes=nsteps,
                            rmax=1,
                            svec=[1.0 2.0 3.0 4.0 5.0],
                            extraargs=QuAPI.QuAPIArgs(),
                            path_integral_routine=QuAPI.build_augmented_propagator)=#

plot(t_TTM.*au2fs, real.(ρ_TTM[:, 1, 1]), label="S0S1")
plot!(t_TTM.*au2fs, real.(ρ_TTM[:, 3, 3]), label="CT")
plot!(t_TTM.*au2fs, real.(ρ_TTM[:, 5, 5]), label="TT")

cond(H0)

eigvals(H0)

fbU = Propagators.calculate_bare_propagators(; Hamiltonian=H0, dt=dt, ntimes=nsteps)
@time t_TTM, ρ_TTM_TEMPO = TTM.propagate(; fbU=fbU,
                            Jw=[Jw],
                            β=β,
                            ρ0=ρ0,
                            dt=dt,
                            ntimes=nsteps,
                            rmax=10,
                            svec=[1.0 2.0 3.0 4.0 5.0],
                            extraargs=TEMPO.TEMPOArgs(),
                            path_integral_routine=TEMPO.build_augmented_propagator)

plot(t_TTM.*au2fs, real.(ρ_TTM_TEMPO[:, 1, 1]), label="S0S1")
plot(t_TTM.*au2fs, real.(ρ_TTM_TEMPO[:, 2, 2]), label="S1S0")
plot!(t_TTM.*au2fs, real.(ρ_TTM_TEMPO[:, 3, 3]), label="CT1")
plot!(t_TTM.*au2fs, real.(ρ_TTM_TEMPO[:, 4, 4]), label="CT2")
plot!(t_TTM.*au2fs, real.(ρ_TTM_TEMPO[:, 5, 5]), label="TT")

plot(tquapi.*au2fs, real.(ρquapi[:, 1, 1]), label="S0S1")
plot!(tquapi.*au2fs, real.(ρquapi[:, 3, 3]), label="CT")
plot!(tquapi.*au2fs, real.(ρquapi[:, 5, 5]), label="TT")

fbU = Propagators.calculate_bare_propagators(; Hamiltonian=H0, dt, ntimes=nsteps)
tq, ρq = TEMPO.propagate(; 
                        fbU=fbU,
                        Jw=[Jw],
                        β=β,
                        ρ0=ρ0,
                        dt,
                        ntimes=nsteps,
                        svec=[1.0 2.0 3.0 4.0 5.0],
                        kmax=3)

plot(tq.*au2fs, real.(ρq[:, 1, 1]), label="S1S0")
plot!(tq.*au2fs, real.(ρq[:, 3, 3]), label="CT")
plot!(tq.*au2fs, real.(ρq[:, 5, 5]), label="TT")

# HEOM calculation here (relies on defining a set of spectral densities (not sure why)

λs = repeat([135.0], 5) * mev2au
γs = repeat([150.0], 5) * mev2au
JwH = Vector{SpectralDensities.DrudeLorentz}()
sys_ops = Vector{Matrix{ComplexF64}}()
for (j, (λ, γ)) in enumerate(zip(λs, γs))
    push!(JwH, SpectralDensities.DrudeLorentz(; λ, γ, Δs=1.0))
    op = zeros(5, 5)
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

plot(times_HEOM.*au2fs, real.(ρs_HEOM[:, 1, 1]), label="S1S0")
plot!(times_HEOM.*au2fs, real.(ρs_HEOM[:, 3, 3]), label="CT")
plot!(times_HEOM.*au2fs, real.(ρs_HEOM[:, 5, 5]), label="TT")

times_BRME, ρs_BRME = BlochRedfield.propagate(;
                                            Hamiltonian=H0,
                                            Jw=JwH,
                                            β,
                                            ρ0,
                                            dt,
                                            ntimes=nsteps,
                                            sys_ops)

plot(times_BRME.*au2fs, real.(ρs_BRME[:, 1, 1]), label="S1S0")
plot!(times_BRME.*au2fs, real.(ρs_BRME[:, 3, 3]), label="CT")
plot!(times_BRME.*au2fs, real.(ρs_BRME[:, 5, 5]), label="TT")


