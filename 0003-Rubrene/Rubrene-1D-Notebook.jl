using QuantumDynamics
using Plots
using LinearAlgebra

const thz2au = 0.0001519828500716
const invcm2au = 4.55633e-6
const au2fs = 0.02418884254
const mev2invcm = 8.066

struct fitsd <: SpectralDensities.AnalyticalSpectralDensity
    ωs :: Vector{Float64}
    jws :: Vector{Float64}
    ωmax :: Real
    Δs :: Real
    classical :: Bool
end

σ = 10*invcm2au
function evaluate(sd::fitsd, ω::Real)
    ωs = sd.ωs
    jws = sd.jws
    
    s = 0.0
    for i in 1:(size(ωs)[1])
        s += jws[i]*exp(-(((ω-ωs[i])/σ)^2))
    end
    s
end

(sd::fitsd)(ω::Real) = evaluate(sd, ω)

ϵ0 = 0.0 # Based on typical NN exciton Hamiltonians
ϵb = 134.0
ϵ2b = -10.7

N = 10

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

ωp = [57.8, 59.6, 89.0, 107.3, 139.1, 639.1, 1011.2, 1344.7, 1593.3] .* invcm2au
ωpg0p = [-1.7, 1.4, 1.6, -0.14, -2.3, -7.5, -3.6, 19.8, -42.0] * invcm2au
g0p = ωpg0p ./ ωp
jws = (g0p.^(2))
wm = maximum(ωp)
Jw = fitsd(ωp, jws, wm, 1.0, false)

ω = 0:0.0001:0.07
plot(ω, Jw.(ω))

nsteps = 5000
ρ0 = Matrix{ComplexF64}(zeros(N, N))
ρ0[5, 5] = 1

T = 50.0:50.0:500.0

β = 1 / (T[6] * 3.16683e-6) # T = 300K

dt = 100 / au2fs

fbU = Propagators.calculate_bare_propagators(; Hamiltonian=H0, dt=dt, ntimes=nsteps)    
t, ρ = TTM.propagate(; fbU=fbU, Jw=[Jw], β=β, ρ0=ρ0, dt=dt, ntimes=nsteps, rmax=1, extraargs=QuAPI.QuAPIArgs(), path_integral_routine=QuAPI.build_augmented_propagator)

## QCPI Run

ω, c = discretize(Jw, 9)
hb = Solvents.HarmonicBath(β, ω, c, [1.0, -1.0], 10)

tc, ρc = QCPI.propagate(; Hamiltonian=H0, Jw, solvent=hb, ρ0, classical_dt=dt / 100, dt, ntimes=nsteps, kmax=1, extraargs=QuAPI.QuAPIArgs(), path_integral_routine=QuAPI.propagate, verbose=false) 

## Running with HEOM

λs = repeat([35.0], N) * invcm2au
γs = 1 ./ (repeat([50.0], N) ./ au2fs)
JwD = Vector{SpectralDensities.DrudeLorentz}()
sys_ops = Vector{Matrix{ComplexF64}}()
for (j, (λ, γ)) in enumerate(zip(λs, γs))
    push!(JwD, SpectralDensities.DrudeLorentz(; λ, γ, Δs=1.0))
    op = zeros(N, N)
    op[j, j] = 1.0
    push!(sys_ops, op)
end

@time th, ρh = HEOM.propagate(; Hamiltonian=H0, ρ0=ρ0, Jw=JwD, β, ntimes=nsteps, dt, sys_ops, num_modes=2, Lmax=3, scaled=true, threshold=1e-10, extraargs=Utilities.DiffEqArgs(; reltol=1e-6, abstol=1e-6))

MSD = []
for i in 1:nsteps
    s = 0.0
    for j in 1:N
        s += real(ρ[i, j, j])*((j-5)^2)
    end
    push!(MSD, s)
end

display("text/plain", MSD)

μ = β * (MSD[nsteps] - MSD[1])/(t[nsteps] - t[1])

t .*= au2fs

plot(t[1:nsteps], MSD, fmt= :png)

plot(t, real.(ρ[:, 5, 5]), label="site 5")
plot!(t, real.(ρ[:, 4, 4]), label="site 4")
plot!(t, real.(ρ[:, 3, 3]), label="site 3")
plot!(t, real.(ρ[:, 2, 2]), label="site 2")
plot!(t, real.(ρ[:, 1, 1]), label="site 1")

plot(th, real.(ρh[:, 5, 5]), label="site 5")
plot!(th, real.(ρh[:, 4, 4]), label="site 4")
plot!(th, real.(ρh[:, 3, 3]), label="site 3")
plot!(th, real.(ρh[:, 2, 2]), label="site 2")
plot!(th, real.(ρh[:, 1, 1]), label="site 1")


