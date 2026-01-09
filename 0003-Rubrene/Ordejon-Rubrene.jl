Threads.nthreads()

using QuantumDynamics
using Plots
using LinearAlgebra

const thz2au = 0.0001519828500716
const invcm2au = 4.55633e-6
const au2fs = 0.02418884254
const mev2invcm = 8.066
const mev2au = mev2invcm * invcm2au
const nm2au = 18.897

ϵ0 = 0.0 # This appears to not change the dynamics at all
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

show(stdout, "text/plain", real.(H0))

# ωp = [57.8, 59.6, 89.0, 107.3, 139.1, 639.1, 1011.2, 1344.7, 1593.3] .* invcm2au
# ωpg0p = [-1.7, 1.4, 1.6, -0.14, -2.3, -7.5, -3.6, 19.8, -42.0] * mev2invcm * invcm2au
# #g0p = ωpg0p ./ ωp
# jws = (ωpg0p .^ 2) ./ ωp
# wm = maximum(ωp)
# Jw = fitsd(ωp, jws, wm, 1.0, false)

# ω = 0:0.0001:0.03
# plot(ω, (Jw.(ω)), xlabel="ω", ylabel="J(ω)", title="Spectral Density", legend=false)

svec = Matrix{Float64}(zeros(1, N))
for i in 1:N
    svec[i] = i
end
svec

Jw = SpectralDensities.DrudeLorentz(λ=135*mev2au, γ=150*mev2au, Δs=1.0)

nsteps = 10000
ρ0 = Matrix{ComplexF64}(zeros(N, N))
ρ0[5, 5] = 1.0

T = 50.0:50.0:500.0

β = 1 / (300 * 3.16683e-6) # T = 300K

dt = 0.1/au2fs # 1fs

fbU = Propagators.calculate_bare_propagators(; Hamiltonian=H0, dt=dt, ntimes=nsteps)
@time t, ρs = TTM.propagate(; fbU=fbU,
                            Jw=[Jw],
                            β=β,
                            ρ0=ρ0,
                            dt=dt,
                            ntimes=nsteps,
                            rmax=1,
                            svec=svec,
                            extraargs=TEMPO.TEMPOArgs(; cutoff=1e-15, maxdim=1000),
                            path_integral_routine=TEMPO.build_augmented_propagator)

for i in 1:N
    plot!(t.*au2fs, real.(ρs[:, i, i]), label="site $i")
end
savefig("rubrene.png")

MSD = []
for i in 1:10000
    s = 0.0
    for j in 1:N
        s += real(ρs[i, j, j])*(j-5)^2
    end
    push!(MSD, s)
end

plot((t[1:10000]).*au2fs, MSD)

# HEOM

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
                                    num_modes=3,
                                    Lmax=5)

plot(times_HEOM.*au2fs, real.(ρs_HEOM[:, 5, 5]), label="site 5")
plot!(times_HEOM.*au2fs, real.(ρs_HEOM[:, 4, 4]), label="site 4")
plot!(times_HEOM.*au2fs, real.(ρs_HEOM[:, 3, 3]), label="site 3")
plot!(times_HEOM.*au2fs, real.(ρs_HEOM[:, 2, 2]), label="site 2")
plot!(times_HEOM.*au2fs, real.(ρs_HEOM[:, 1, 1]), label="site 1")

MSD_HEOM = []
for i in 1:30000
    s = 0.0
    for j in 1:N
        s += real(ρs_HEOM[i, j, j])*(j-5)^2
    end
    push!(MSD_HEOM, s)
end

plot((times_HEOM[1:30000]).*au2fs, MSD_HEOM)

μ = (MSD_HEOM[30000] - MSD_HEOM[1])/((times_HEOM[30000] - times_HEOM[1])*au2fs)


