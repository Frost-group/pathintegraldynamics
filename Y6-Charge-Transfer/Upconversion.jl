# Simple script to execute HEOM and TEMPO-TTM for molecular trimer Singlet Fission simulations.

using QuantumDynamics
using LinearAlgebra
using DelimitedFiles
using OrdinaryDiffEq


include("hamiltonian-models.jl")
include("ratematrix-solver.jl")

const thz2au = 0.0001519828500716
const invcm2au = 4.55633e-6
const au2fs = 0.02418884254
const mev2invcm = 8.066
const mev2au = mev2invcm * invcm2au
const nm2au = 18.897


#= For now : use code upconversion-ratematrix.jl

function UpconvertRateMatrix()
    λs, γ, H0 = Y6UpconversionDimerHamiltonianWTriplet(162.0, 71.0)

    N = 5

    β = 1052.0 # 300K
    ts, ps = propagate(N, H0, λs, β, 4000, 0.25/au2fs)
    
    ts = ts.*au2fs
    
   # println(size(ts))
    #println(size(ps[:, 1]))
    #println(ps)
    #println(ts)

    open("upconversion-ratematrix-populations.stdout", "w") do io
       pops = [real.(ps[:, i]) for i in 1:N]
       tpops = [ts pops...]
       writedlm(io, tpops, ' ')
   end
end
=#

function Boltzmann()
    
    λ, γ, H = Y6UpconversionDimerHamiltonianWTriplet(71.0, 162.0)
    
    N = length(λ)

    E = [real(H[i, i]) for i in 1:N]
    
    β = 1052.0

    Z = sum(ℯ.^(-β*E))

    pops = [ℯ^(-β*Ei)/Z for Ei in E]
    

    open("boltzmann-populations.stdout", "w") do io
        writedlm(io, pops, ' ')
    end

end


"""
UpconversionHEOM

"""


function UpconversionHEOM(s, V; dt=0.25/au2fs, nsteps=4000, L=4, K=2)

    λs, γs, H0 = Y6UpconversionDimerHamiltonian(s, s, V, V)

    N = length(λs)

    ρ0 = Matrix{ComplexF64}(zeros(N, N))
    ρ0[1, 1] = 1.0
    β = 1 / (300 * 3.16683e-6) # T = 300K
    svec = Matrix{Float64}(zeros(1, N))
    for i in 1:N
        svec[i] = i
    end


    JwH = Vector{SpectralDensities.DrudeLorentz}()
    sys_ops = Vector{Matrix{ComplexF64}}()
    for (j, (λ, γ)) in enumerate(zip(λs, γs))
       push!(JwH, SpectralDensities.DrudeLorentz(; λ, γ, Δs=1.0))
       op = zeros(N, N)
       op[j, j] = 1.0
       push!(sys_ops, op)
    end

    times_HEOM, ρs = HEOM.propagate(;
                                    Hamiltonian=H0,
                                    ρ0,
                                    β,
                                    dt,
                                    ntimes=nsteps,
                                    Jw=JwH,
                                    sys_ops=sys_ops,
                                    num_modes=K,
                                    Lmax=L)


    open("upconversion-dimer-populations-heom.stdout", "w") do io
       pops = [real.(ρs[:, i, i]) for i in 1:N]
       tpops = [times_HEOM pops...]
       writedlm(io, tpops, ' ')
   end

end


"""
UpconversionRedfield

"""


function UpconversionRedfield(s, V; dt=0.25/au2fs, nsteps=4000)
   
    λs, γs, H0 = Y6UpconversionDimerHamiltonian(s, s, V, V)


    N = length(λs)


    ρ0 = Matrix{ComplexF64}(zeros(N, N))
    ρ0[1, 1] = 1.0
    β = 1 / (300 * 3.16683e-6) # T = 300K

    λ_av = sum(λs)/length(λs)
    γ_av = sum(γs)/length(γs)

    η = max(2*λ_av/(β*γ_av^2), 2*λ_av/(π*γ_av))
    
    println(η)

    if η > 1
        @info "η = $η > 1 ; Redfield is inaccurate here"
    end
    svec = Matrix{Float64}(zeros(1, N))
    for i in 1:N
        svec[i] = i
    end


    JwH = Vector{SpectralDensities.DrudeLorentz}()
    sys_ops = Vector{Matrix{ComplexF64}}()
    for (j, (λ, γ)) in enumerate(zip(λs, γs))
       push!(JwH, SpectralDensities.DrudeLorentz(; λ, γ, Δs=1.0))
       op = zeros(N, N)
       op[j, j] = 1.0
       push!(sys_ops, op)
    end

    times_BRME, ρs = BlochRedfield.propagate(;
                                    Hamiltonian=H0,
                                    ρ0,
                                    β,
                                    dt,
                                    ntimes=nsteps,
                                    Jw=JwH,
                                    sys_ops=sys_ops,
                                   )

   # open("stdout-files4/upconversion-populations-brme-s-$s-V-$V-t-$E.stdout", "w") do io
   open("upconversion-redfield-dimer.stdout", "w") do io
       pops = [real.(ρs[:, i, i]) for i in 1:N]
       tpops = [times_BRME pops...]
       writedlm(io, tpops, ' ')
   end

end


"""
UpconversionTTM

"""

function UpconversionTTM(; dt=0.25/au2fs, nsteps=4000, rmax=15)
    

    reorg = 157 * mev2au

    cutoff = 1600 * invcm2au

    Efe = 2000.0
    Ect(r) = (2.19 -  4.959/(r))*1000  # Enter r in angstrom ; best fit equation for Ect
    Dh = 55.7
    De = 72.0
    V = -76.0
    
    Ec = Ect(9.29)
    
    N = 4



    H0 = Matrix{ComplexF64}([    ## Check this pls
        Efe V Dh De
        V Efe De Dh
        Dh De Ec 0.0
        De Dh 0.0 Ec
    ]) * mev2au

    ρ0 = Matrix{ComplexF64}(zeros(N, N))
    ρ0[2, 2] = 1.0
    β = 1 / (300 * 3.16683e-6) # T = 300K
    
    svec = [1.0 2.0 3.0 4.0]
    
    

    Jw = SpectralDensities.DrudeLorentz(λ=reorg, γ=cutoff, Δs=1.0)

    fbU = Propagators.calculate_bare_propagators(; Hamiltonian=H0, dt=dt, ntimes=nsteps)

    ts, ρs = TTM.propagate(; fbU=fbU,
                            Jw=[Jw],
                            β=β,
                            ρ0=ρ0,
                            dt=dt,
                            ntimes=nsteps,
                            rmax=rmax,
                            svec=svec,
                            extraargs=TEMPO.TEMPOArgs(),
                            path_integral_routine=TEMPO.build_augmented_propagator)
    
    open("upconversion-populations-ttm.stdout", "w") do io
        pops = [real.(ρs[:, i, i]) for i in 1:N]
        tpops = [ts pops...]
        writedlm(io, tpops, ' ')
    end
end


#UpconvertRateMatrix()

#UpconversionRedfield(71.0, 162.0)

UpconversionHEOM(71.0, 162.0)

