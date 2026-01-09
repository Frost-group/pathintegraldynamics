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
v = 50*invcm2au

N = 8

H0 = Matrix{ComplexF64}(zeros(N, N))

for i in 1:N
    H0[i,i] = ϵ0
    if i <= N-1
        H0[i, i+1] = v
    end
    if i >= 2
        H0[i, i-1] = v
    end
end

show(stdout, "text/plain", real.(H0))

nsteps = 40000
ρ0 = Matrix{ComplexF64}(zeros(N, N))
ρ0[1, 1] = 1.0

T = 50.0:50.0:500.0

β = 1 / (300 * 3.16683e-6) # T = 300K

dt = 0.25/au2fs

svec = Matrix{Float64}(zeros(1, N))
for i in 1:N
    svec[i] = i
end
svec

Jw = SpectralDensities.DrudeLorentz(λ=161.5*invcm2au, γ=41*invcm2au, Δs=1.0)

ω = 0.0000001:0.0001:0.07
plot(ω, Jw.(ω))

using DelimitedFiles

open("Jw_1", "w") do io
    writedlm(io, [ω  Jw.(ω)], ' ')
end

fbU = Propagators.calculate_bare_propagators(; Hamiltonian=H0, dt=dt, ntimes=nsteps)
@time t, ρs = TTM.propagate(; fbU=fbU,
                            Jw=[Jw],
                            β=β,
                            ρ0=ρ0,
                            dt=dt,
                            ntimes=nsteps,
                            rmax=7,
                            svec=svec,
                            extraargs=TEMPO.TEMPOArgs(; cutoff=1e-15, maxdim=1000),
                            path_integral_routine=TEMPO.build_augmented_propagator)

plot(t.*au2fs, real.(ρs[:, 1, 1]), label="site 1")
plot!(t.*au2fs, real.(ρs[:, 2, 2]), label="site 2")
plot!(t.*au2fs, real.(ρs[:, 3, 3]), label="site 3")
plot!(t.*au2fs, real.(ρs[:, 4, 4]), label="site 4")
plot!(t.*au2fs, real.(ρs[:, 5, 5]), label="site 5")
plot!(t.*au2fs, real.(ρs[:, 6, 6]), label="site 6")
plot!(t.*au2fs, real.(ρs[:, 7, 7]), label="site 7")
plot!(t.*au2fs, real.(ρs[:, 8, 8]), label="site 8")

MSD = []
for i in 1:40000
    s = 0.0
    for j in 1:N
        s += real(ρs[i, j, j])*(j-1)^2
    end
    push!(MSD, s)
end

plot((t[2:40000]).*au2fs, MSD[2:40000], xscale=:log10)

dMSD_dt = [(MSD[i] - MSD[i-1])/(dt*au2fs) for i in 2:40000]

plot((t[2:40000]).*au2fs, dMSD_dt, xscale=:log10)

λs = repeat([161.5], N) * invcm2au
γs = repeat([41.0], N) * invcm2au
JwH = Vector{SpectralDensities.DrudeLorentz}()
sys_ops = Vector{Matrix{ComplexF64}}()
for (j, (λ, γ)) in enumerate(zip(λs, γs))
    push!(JwH, SpectralDensities.DrudeLorentz(; λ, γ, Δs=1.0))
    op = zeros(N, N)
    op[j, j] = 1.0
    push!(sys_ops, op)
end


@time times_HEOM, ρs_HEOM = HEOM.propagate(;
                                    Hamiltonian=H0,
                                    ρ0,
                                    β,
                                    dt,
                                    ntimes=nsteps,
                                    Jw=JwH,
                                    sys_ops=sys_ops,
                                    num_modes=1,
                                    Lmax=5)

plot(times_HEOM.*au2fs, real.(ρs_HEOM[:, 1, 1]), label="site 1", xscale=:log10)
plot!(times_HEOM.*au2fs, real.(ρs_HEOM[:, 2, 2]), label="site 2", xscale=:log10)
plot!(times_HEOM.*au2fs, real.(ρs_HEOM[:, 3, 3]), label="site 3", xscale=:log10)
plot!(times_HEOM.*au2fs, real.(ρs_HEOM[:, 4, 4]), label="site 4", xscale=:log10)
plot!(times_HEOM.*au2fs, real.(ρs_HEOM[:, 5, 5]), label="site 5", xscale=:log10)
plot!(times_HEOM.*au2fs, real.(ρs_HEOM[:, 6, 6]), label="site 6", xscale=:log10)
plot!(times_HEOM.*au2fs, real.(ρs_HEOM[:, 7, 7]), label="site 7", xscale=:log10)
plot!(times_HEOM.*au2fs, real.(ρs_HEOM[:, 8, 8]), label="site 8", xscale=:log10)

MSD_HEOM = []
for i in 1:10000
    s = 0.0
    for j in 1:N
        s += real(ρs_HEOM[i, j, j])*(j-1)^2
    end
    push!(MSD_HEOM, s)
end

plot((times_HEOM[1:40000]).*au2fs, MSD_HEOM)

dMSD_dt = [(MSD_HEOM[i] - MSD_HEOM[i-1])/(dt*au2fs) for i in 2:40000]

plot((times_HEOM[1:40000]).*au2fs, dMSD_dt)

using DelimitedFiles

dlm = readdlm("AMC-holstein-output-rmax4", Float64)

t = dlm[:, 1]
t.*=au2fs
p1 = dlm[:, 2]
p2 = dlm[:, 5]
p3 = dlm[:, 10]
p4 = dlm[:, 17]
p5 = dlm[:, 26]
p6 = dlm[:, 37]
p7 = dlm[:, 50]
p8 = dlm[:, 65]

plot(t, p1, label="site 1",title="SMatPI", xscale=:log10)
plot!(t, p2, label="site 2",title="SMatPI")
plot!(t, p3, label="site 3",title="SMatPI")
plot!(t, p4, label="site 4",title="SMatPI")
#plot!(t, p5, label="site 5",title="SMatPI")
#plot!(t, p6, label="site 6",title="SMatPI")
#plot!(t, p7, label="site 7",title="SMatPI")
#plot!(t, p8, label="site 8",title="SMatPI")

MSD_smatpi = []
for i in 1:400000
    s = (p1[i]*0 + p2[i]*1 + p3[i]*4 + p4[i]*9 + p5[i]*16 + p6[i]*25 + p7[i]*36 + p8[i]*49)*(nm2au)^2*0.25
    push!(MSD_smatpi, s)
end

plot(t, MSD_smatpi, xscale=:log10)

dMSD_dt = [(MSD_smatpi[i] - MSD_smatpi[i-1])/(dt*au2fs) for i in 2:400000]

plot((t[2:400000]), dMSD_dt, xscale=:log10)


