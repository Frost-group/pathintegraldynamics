using QuantumDynamics
using Plots
using LinearAlgebra

const thz2au = 0.0001519828500716
const invcm2au = 4.55633e-6
const au2fs = 0.02418884254
const mev2invcm = 8.066
const mev2au = mev2invcm * invcm2au
const nm2au = 18.897

N=2

H0 = Matrix{ComplexF64}([
        0.0 50.0
        50.0 -75.0
    ]) * mev2au

λs = repeat([25.0], N) * mev2au
γs = repeat([150.0], N) * mev2au
JwH = Vector{SpectralDensities.DrudeLorentz}()
sys_ops = Vector{Matrix{ComplexF64}}()
for (j, (λ, γ)) in enumerate(zip(λs, γs))
    push!(JwH, SpectralDensities.DrudeLorentz(; λ, γ, Δs=1.0))
    op = zeros(N, N)
    op[j, j] = 1.0
    push!(sys_ops, op)
end

dt = 1/au2fs 
nsteps = 150
ρ0 = Matrix{ComplexF64}(zeros(N, N))
ρ0[1, 1] = 1.0

β = 1 / (300 * 3.16683e-6) # T = 300K

@time t_HEOM, ρs_HEOM = HEOM.propagate(;
                                    Hamiltonian=H0,
                                    ρ0,
                                    β,
                                    dt,
                                    ntimes=nsteps,
                                    Jw=JwH,
                                    sys_ops=sys_ops,
                                    num_modes=2,
                                    Lmax=3)

Jw = SpectralDensities.DrudeLorentz(λ=25*mev2au, γ=150*mev2au, Δs=1.0)
fbU = Propagators.calculate_bare_propagators(; Hamiltonian=H0, dt=dt, ntimes=nsteps)
@time t_TTM, ρ_TTM = TTM.propagate(; fbU=fbU,
                            Jw=[Jw],
                            β=β,
                            ρ0=ρ0,
                            dt=dt,
                            ntimes=nsteps,
                            rmax=10,
                            svec=[1.0 2.0],
                            extraargs=TEMPO.TEMPOArgs(),
                            path_integral_routine=TEMPO.build_augmented_propagator)

@time tq, ρq = QuAPI.propagate(; 
                        fbU=fbU,
                        Jw=[Jw],
                        β=β,
                        ρ0=ρ0,
                        dt,
                        ntimes=nsteps,
                        svec=[1.0 2.0],
                        kmax=7)

@time t_BRME, ρs_BRME = BlochRedfield.propagate(;
                                            Hamiltonian=H0,
                                            Jw=JwH,
                                            β,
                                            ρ0,
                                            dt,
                                            ntimes=nsteps,
                                            sys_ops)

plot(t_HEOM.*au2fs, real.(ρs_HEOM[:,1,1]), label="HEOM")
plot!(t_TTM.*au2fs, real.(ρ_TTM[:,1,1]), label="TTM")
plot!(t_BRME.*au2fs, real.(ρs_BRME[:,1,1]), label="BRME")
plot!(tq.*au2fs, real.(ρq[:,1,1]), label="QuAPI")

ES = 250.0
ECT = 200.0
ETT = 0.0

VSCT = 3.0
VCTTT = 17

N = 3

H_seq = Matrix{ComplexF64}([
        ES VSCT 0.0
        VSCT ECT VCTTT
        0.0 VCTTT ETT
    ])*mev2au

λs = repeat([25.0], N) * mev2au
γs = repeat([150.0], N) * mev2au
JwH = Vector{SpectralDensities.DrudeLorentz}()
sys_ops = Vector{Matrix{ComplexF64}}()
for (j, (λ, γ)) in enumerate(zip(λs, γs))
    push!(JwH, SpectralDensities.DrudeLorentz(; λ, γ, Δs=1.0))
    op = zeros(N, N)
    op[j, j] = 1.0
    push!(sys_ops, op)
end

dt = 1/au2fs 
nsteps = 14000
ρ0 = Matrix{ComplexF64}(zeros(N, N))
ρ0[1, 1] = 1.0

β = 1 / (300 * 3.16683e-6) # T = 300K

@time t_HEOM, ρs_HEOM = HEOM.propagate(;
                                    Hamiltonian=H_seq,
                                    ρ0,
                                    β,
                                    dt,
                                    ntimes=nsteps,
                                    Jw=JwH,
                                    sys_ops=sys_ops,
                                    num_modes=2,
                                    Lmax=3)

#=Jw = SpectralDensities.DrudeLorentz(λ=25*mev2au, γ=150*mev2au, Δs=1.0)
fbU = Propagators.calculate_bare_propagators(; Hamiltonian=H_seq, dt=dt, ntimes=nsteps)
@time t_TTM, ρ_TTM = TTM.propagate(; fbU=fbU,
                            Jw=[Jw],
                            β=β,
                            ρ0=ρ0,
                            dt=dt,
                            ntimes=nsteps,
                            rmax=10,
                            svec=[1.0 2.0 3.0],
                            extraargs=TEMPO.TEMPOArgs(),
                            path_integral_routine=TEMPO.build_augmented_propagator)=#

Jw = SpectralDensities.DrudeLorentz(; λ=25*mev2au, γ=150*mev2au, Δs=1.0)
ωs, cs = SpectralDensities.discretize(Jw, 20)
hb = Solvents.HarmonicBath(β, ωs, cs, [1.0, 2.0, 3.0], 1000)
@time tc, ρc = QCPI.propagate(; Hamiltonian=H0, Jw=Jw, solvent=hb, ρ0, classical_dt=dt / 100, dt, ntimes=nsteps, kmax=3, svec=[1.0 2.0 3.0], extraargs=QuAPI.QuAPIArgs(), path_integral_routine=QuAPI.propagate, verbose=true)

@time t_BRME, ρs_BRME = BlochRedfield.propagate(;
                                            Hamiltonian=H_seq,
                                            Jw=JwH,
                                            β,
                                            ρ0,
                                            dt,
                                            ntimes=nsteps,
                                            sys_ops)

plot(t_HEOM.*au2fs, real.(ρs_HEOM[:,1,1]), label="HEOM S1", ls=:dot, lc=:black)
plot!(t_HEOM.*au2fs, real.(ρs_HEOM[:,2,2]), label="HEOM CT", ls=:dot, lc=:black)
plot!(t_HEOM.*au2fs, real.(ρs_HEOM[:,3,3]), label="HEOM TT", ls=:dot, lc=:black)
#plot!(t_TTM.*au2fs, real.(ρ_TTM[:,1,1]), label="TTM", lc=:blue)
#plot!(t_TTM.*au2fs, real.(ρ_TTM[:,2,2]), label="TTM", lc=:blue)
#plot!(t_TTM.*au2fs, real.(ρ_TTM[:,3,3]), label="TTM", lc=:blue)
plot!(t_BRME.*au2fs, real.(ρs_BRME[:,1,1]), label="BRME", lc=:green)
plot!(t_BRME.*au2fs, real.(ρs_BRME[:,2,2]), label="BRME", lc=:green)
plot!(t_BRME.*au2fs, real.(ρs_BRME[:,3,3]), label="BRME", lc=:green)

ES = 80.0
ECT = 330.0
ETT = 0.0

VSCT = 30.0
VCTTT = 30.0

N = 3

H_sup = Matrix{ComplexF64}([
        ES VSCT 0.0
        VSCT ECT VCTTT
        0.0 VCTTT ETT
    ])*mev2au

λs = repeat([25.0], N) * mev2au
γs = repeat([150.0], N) * mev2au
JwH = Vector{SpectralDensities.DrudeLorentz}()
sys_ops = Vector{Matrix{ComplexF64}}()
for (j, (λ, γ)) in enumerate(zip(λs, γs))
    push!(JwH, SpectralDensities.DrudeLorentz(; λ, γ, Δs=1.0))
    op = zeros(N, N)
    op[j, j] = 1.0
    push!(sys_ops, op)
end

dt = 1/au2fs 
nsteps = 1400
ρ0 = Matrix{ComplexF64}(zeros(N, N))
ρ0[1, 1] = 1.0

β = 1 / (300 * 3.16683e-6) # T = 300K

@time t_HEOM, ρs_HEOM = HEOM.propagate(;
                                    Hamiltonian=H_sup,
                                    ρ0,
                                    β,
                                    dt,
                                    ntimes=nsteps,
                                    Jw=JwH,
                                    sys_ops=sys_ops,
                                    num_modes=3,
                                    Lmax=5)

#=Jw = SpectralDensities.DrudeLorentz(λ=25*mev2au, γ=150*mev2au, Δs=1.0)
fbU = Propagators.calculate_bare_propagators(; Hamiltonian=H_seq, dt=dt, ntimes=nsteps)
@time t_TTM, ρ_TTM = TTM.propagate(; fbU=fbU,
                            Jw=[Jw],
                            β=β,
                            ρ0=ρ0,
                            dt=dt,
                            ntimes=nsteps,
                            rmax=10,
                            svec=[1.0 2.0 3.0],
                            extraargs=TEMPO.TEMPOArgs(),
                            path_integral_routine=TEMPO.build_augmented_propagator)=#

@time t_BRME, ρs_BRME = BlochRedfield.propagate(;
                                            Hamiltonian=H_sup,
                                            Jw=JwH,
                                            β,
                                            ρ0,
                                            dt,
                                            ntimes=nsteps,
                                            sys_ops)

# plot(t_HEOM.*au2fs, real.(ρs_HEOM[:,1,1]), label="HEOM S1", ls=:dot, lc=:black)
# plot!(t_HEOM.*au2fs, real.(ρs_HEOM[:,2,2]), label="HEOM CT", ls=:dot, lc=:black)
# plot!(t_HEOM.*au2fs, real.(ρs_HEOM[:,3,3]), label="HEOM TT", ls=:dot, lc=:black)
# #plot!(t_TTM.*au2fs, real.(ρ_TTM[:,1,1]), label="TTM", lc=:blue)
# #plot!(t_TTM.*au2fs, real.(ρ_TTM[:,2,2]), label="TTM", lc=:blue)
# #plot!(t_TTM.*au2fs, real.(ρ_TTM[:,3,3]), label="TTM", lc=:blue)
# plot!(t_BRME.*au2fs, real.(ρs_BRME[:,1,1]), label="BRME", lc=:green)
# plot!(t_BRME.*au2fs, real.(ρs_BRME[:,2,2]), label="BRME", lc=:green)
# plot!(t_BRME.*au2fs, real.(ρs_BRME[:,3,3]), label="BRME", lc=:green)

using Profile, FlameGraphs, ImageShow

ET = 1720.0
ECT = 250.0 + ET
ES = 500.0 + ET


# Sequential

Hpd = Matrix{ComplexF64}([
        ECT 116.0 -152.0 145.0 0.0
        116.0  ES   0.0 0.0 145.0
        -152.0 0.0 ET 0 133.0
        145.0 0.0 0 ES 116.0
        0.0 145.0 133.0 116.0 ECT
    ]) * mev2au

N = 5

λs = repeat([25.0], N) * mev2au
γs = repeat([150.0], N) * mev2au
JwH = Vector{SpectralDensities.DrudeLorentz}()
sys_ops = Vector{Matrix{ComplexF64}}()
for (j, (λ, γ)) in enumerate(zip(λs, γs))
    push!(JwH, SpectralDensities.DrudeLorentz(; λ, γ, Δs=1.0))
    op = zeros(N, N)
    op[j, j] = 1.0
    push!(sys_ops, op)
end

dt = 1/au2fs 
nsteps = 500
ρ0 = Matrix{ComplexF64}(zeros(N, N))
ρ0[2, 2] = 1.0

β = 1 / (300 * 3.16683e-6) # T = 300K

@time t_HEOM, ρs_HEOM = HEOM.propagate(;
                                    Hamiltonian=Hpd,
                                    ρ0,
                                    β,
                                    dt,
                                    ntimes=nsteps,
                                    Jw=JwH,
                                    sys_ops=sys_ops,
                                    num_modes=3,
                                    Lmax=5)

plot(t_HEOM.*au2fs, real.(ρs_HEOM[:,1,1]), label="HEOM CT", lc=:orange)
plot!(t_HEOM.*au2fs, real.(ρs_HEOM[:,2,2]), label="HEOM S1", lc=:blue)
plot!(t_HEOM.*au2fs, real.(ρs_HEOM[:,3,3]), label="HEOM TT", lc=:green)

Profile.clear()

@profile t_BRME, ρs_BRME = BlochRedfield.propagate(;
                                            Hamiltonian=Hpd,
                                            Jw=JwH,
                                            β,
                                            ρ0,
                                            dt,
                                            ntimes=nsteps,
                                            sys_ops)

g = flamegraph(C=true)

img = flamepixels(g)

plot(t_BRME.*au2fs, real.(ρs_BRME[:,1,1]), label="BRME CT", lc=:orange)
plot!(t_BRME.*au2fs, real.(ρs_BRME[:,2,2]), label="BRME S1", lc=:blue)
plot!(t_BRME.*au2fs, real.(ρs_BRME[:,3,3]), label="BRME TT", lc=:green)

Jw = SpectralDensities.DrudeLorentz(λ=25*mev2au, γ=150*mev2au, Δs=1.0)
fbU = Propagators.calculate_bare_propagators(; Hamiltonian=Hpd, dt=dt, ntimes=nsteps)
@time t_TTM, ρ_TTM = TTM.propagate(; fbU=fbU,
                            Jw=[Jw],
                            β=β,
                            ρ0=ρ0,
                            dt=dt,
                            ntimes=nsteps,
                            rmax=15,
                            svec=[1.0 2.0 3.0 4.0 5.0],
                            extraargs=TEMPO.TEMPOArgs(),
                            path_integral_routine=TEMPO.build_augmented_propagator)

plot(t_TTM.*au2fs, real.(ρ_TTM[:,1,1]), label="TTM CT", lc=:orange)
plot!(t_TTM.*au2fs, real.(ρ_TTM[:,2,2]), label="TTM S1", lc=:blue)
plot!(t_TTM.*au2fs, real.(ρ_TTM[:,3,3]), label="TTM TT", lc=:green)

using DelimitedFiles
ω = 0:0.0001:0.5
writedlm("Jw_1", zip(ω, Jw.(ω)))

show(stdout, "text/plain", real.(Hpd))

eigvals(Hpd)

cond(Hpd)
