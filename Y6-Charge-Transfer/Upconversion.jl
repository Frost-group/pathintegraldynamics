# Simple script to execute HEOM and TEMPO-TTM for molecular trimer Singlet Fission simulations.

#using QuantumDynamics
using LinearAlgebra
using DelimitedFiles
using OrdinaryDiffEq

include("ratematrix-solver.jl")

const thz2au = 0.0001519828500716
const invcm2au = 4.55633e-6
const au2fs = 0.02418884254
const mev2invcm = 8.066
const mev2au = mev2invcm * invcm2au
const nm2au = 18.897


"""
Y6UpconversionTrimerHamiltonian

"""


function Y6UpconversionTrimerHamiltonian()
    
    # D1, D3, D8 trimer dynamics - FE vs CT character over time.

    λ = [repeat([157*mev2au], 3)..., repeat([240*mev2au], 6)...]

    γ = repeat([1600 * invcm2au], 9)

    Ef1 = 2046.0 # eV
    Ef3 = 2046.0
    Ef8 = 2046.0

    Ect(r) = (2.19 -  4.959/(r))*1000  # Enter r in angstrom ; best fit equation for Ect
   
    r1 = 9.29
    r3 = 13.84
    r8 = 15.97
    Dh1 = 55.7
    Dh3 = -27.3
    Dh8 = 0.0
    De1 = 72.0
    De3 = -15.0
    De8 = 0.0
   
    V1 = -76.0
    V3 = -6.1
    V8 = -9.0
   
    Ec1 = Ect(r1) # 1656.2 eV
    Ec3 = Ect(r3) # 1831.69 eV
    Ec8 = Ect(r8) # 1879.48 eV
   
    te1 = 66.1
    te3 = -27.0
    te8 = 0.0
    th1 = 49.6
    th3 = -11.3
    th8 = 0.0

    
    H0 = Matrix{ComplexF64}([
         Ef1 V1 V8 Dh1 Dh8 De1 0.0 De8 0.0
         V1 Ef3 V3 De1 0.0 Dh1 Dh3 0.0 De3
         V8 V3 Ef8 0.0 De8 0.0 De3 Dh8 Dh3
         Dh1 De1 0.0 Ec1 th3 0.0 0.0 0.0 te8
         Dh8 0.0 De8 th3 Ec8 0.0 te1 0.0 0.0  
         De1 Dh1 0.0 0.0 0.0 Ec1 th8 te3 0.0
         0.0 Dh3 De3 0.0 te1 th8 Ec3 0.0 0.0
         De8 0.0 Dh8 0.0 0.0 te3 0.0 Ec8 th1
         0.0 De3 Dh3 te8 0.0 0.0 0.0 th1 Ec3
        ])
    
    return λ, γ, H0

end

function Y6UpconversionDimerHamiltonianWTriplet(soc::Float64, Vct::Float64)
    
    # Params from Samuele paper

    # D1 dimer dynamics - FE vs CT character over time.
    
    reorg = [157.0, 157.0, 240.0, 240.0, 157.0]*mev2au

    cutoff = repeat([1600 * invcm2au], 5)

    #Efe = 2046.0
    Ect(r) = (2.19 -  4.959/(r))*1000  # Enter r in angstrom ; best fit equation for Ect
    Dh = 55.7
    De = 72.0
    V = -76.0
   
    Ec = Ect(9.29)
   
    # Dimer Hamiltonian with singlet and CT states

    #=H0 = Matrix{ComplexF64}([
        Efe V Dh De
        V Efe De Dh
        Dh De Ec 0.0
        De Dh 0.0 Ec
    ]) * mev2au =#

    
    # Zhenghan triplet values
#    Et1 = 1350.0
#    Et2 = 1393.0
    
#    Ett = Et1 + Et2


    # Dimer Hamiltonian including triplet with placeholder SoC
    
    # Zhenghan singlet values

    Efe1 = 1872.0
    Efe2 = 1886.0


    Ett = 1350.0

    H0 = Matrix{ComplexF64}([
        Efe1 V Dh De soc 
        V Efe2 De Dh soc 
        Dh De Ec 0.0 Vct 
        De Dh 0.0 Ec Vct 
        soc soc Vct Vct Ett
    ]) * mev2au 


    return reorg, cutoff, H0
end

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


function UpconversionHEOM(; dt=0.25/au2fs, nsteps=4000, L=1, K=2)

    soc = 71.0
    Vct = 162.0

    #λs, γs, H0 = Y6UpconversionDimerHamiltonianWTriplet(soc, Vct)
    λs, γs, H0 = Y6UpconversionTrimerHamiltonian()

    N = length(λs)

    ρ0 = Matrix{ComplexF64}(zeros(N, N))
    ρ0[2, 2] = 1.0
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


    open("upconversion-trimer-populations-heom.stdout", "w") do io
       pops = [real.(ρs[:, i, i]) for i in 1:N]
       tpops = [times_HEOM pops...]
       writedlm(io, tpops, ' ')
   end

end


"""
UpconversionRedfield

"""


function UpconversionRedfield(s, V; dt=0.25/au2fs, nsteps=4000)
   
    #λs, γs, H0 = Y6UpconversionDimerHamiltonianWTriplet(s, V)

    λs, γs, H0 = Y6UpconversionTrimerHamiltonian()
    N = length(λs)

    ρ0 = Matrix{ComplexF64}(zeros(N, N))
    ρ0[2, 2] = 1.0
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
                                    extraargs = Utilities.DiffEqArgs(reltol=1e-10, abstol=1e-10, solver=RadauIIA3)
                                   )

   # open("stdout-files4/upconversion-populations-brme-s-$s-V-$V-t-$E.stdout", "w") do io
   open("upconversion-redfield-trimer.stdout", "w") do io
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


UpconvertRateMatrix()

# UpconversionRedfield(s, V, E) - inputs - SOC (mev), CT-TT coupling(mev), Triplet pair energy(mev)

#UpconversionRedfield(0.03, 0.1, 2000.0)

#UpconversionRedfield(71.0, 162.0)

#=for i in 1:5:100
    for j in 2:8:200
        UpconversionRedfield(i*1.0, j*1.0, 2743.0)
    end
end=#

#UpconversionTTM() 

#UpconversionHEOM() 


