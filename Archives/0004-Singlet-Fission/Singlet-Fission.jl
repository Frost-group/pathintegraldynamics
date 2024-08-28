# Five State Singlet Fission Model with Redfield, HEOM and QuAPI/TEMPO-TTM.


using QuantumDynamics
using DelimitedFiles

const thz2au = 0.0001519828500716
const invcm2au = 4.55633e-6
const au2fs = 0.02418884254
const mev2invcm = 8.066
const mev2au = mev2invcm * invcm2au
const nm2au = 18.897


function SingletFissionDimerModel(; dt=0.5/au2fs, nsteps=200)
# Grozema PDI Dimer

#   Es = 2.22 * 1000
#   Ec = 3.11 * 1000
#   Et =  1.05 * 1000
#   t2e = -0.42
#   thh = -125.0
#   tll = 145.0
#   thl = -125.0
#   
#   H0 = Matrix{ComplexF64}([
#       Es 0 tll -thh t2e
#       0 Es -thh tll t2e
#       tll -thh Ec 0 sqrt(3/2)*thl
#       -thh tll 0 Ec sqrt(3/2)*thl
#       t2e t2e sqrt(3/2)*thl sqrt(3/2)*thl 2*Et
#                           ]) * mev2au

#   λ = 0.135 * 1000 * mev2au
#   γ = 0.15 * 1000 * mev2au
    
# Troisi Rubrene Dimer

    H0 = Matrix{ComplexF64}([
        4.018 0.058 -0.086 0.168 0.000
        0.058 4.017 -0.168 0.086 0.000
        -0.086 -0.168 5.043 0.000 0.000
        0.168 0.086 0.000 5.042 0.000
        0.000 0.000 0.000 0.000 4.006
                             ]) * 1000 * mev2au

    λ = 330*mev2au
    γ = 1500*invcm2au
    return H0, λ, γ
end


function SingletFissionRedfield(; dt=0.5/au2fs, nsteps=200)
    H0, λ, γ = SingletFissionDimerModel()

    N = 5 
    
    λs = repeat([λ], N)
    γs = repeat([γ], N)

    ρ0 = Matrix{ComplexF64}(zeros(N, N))
    ρ0[1, 1] = 1.0

    T = 50.0:50.0:500.0

    β = 1 / (300 * 3.16683e-6) # T = 300K

    η = max(2*λ/(β*γ^2), 2*λ/(π*γ))
    
    println(η)
    if η > 1
        @warn "η = $η > 1 , Redfield may not be accurate"
    end
    
    

    Jw = Vector{SpectralDensities.DrudeLorentz}()
    sys_ops = Vector{Matrix{ComplexF64}}()
    for (j, (λ, γ)) in enumerate(zip(λs, γs))
        push!(Jw, SpectralDensities.DrudeLorentz(; λ, γ, Δs=1.0))
        op = zeros(N, N)
        op[j, j] = 1.0
        push!(sys_ops, op)
    end

    ts, ρs = BlochRedfield.propagate(;
                               Hamiltonian=H0,
                               Jw=Jw,
                               β,
                               ρ0,
                               dt,
                               ntimes=nsteps,
                               sys_ops)
    
    open("redfield-SF-populations-Troisi-Rubrene.txt", "w") do io
        pops = [real.(ρs[:, i, i]) for i in 1:N]
        tpops = [ts pops...]
        writedlm(io, tpops, ' ')
    end


end


function SingletFissionHEOM(; dt=0.5/au2fs, nsteps=200, Lmax=5, num_modes=2)
    H0, λ, γ = SingletFissionDimerModel()
    
    N = 5

    
    λs = repeat([λ], N)
    γs = repeat([γ], N)
    ρ0 = Matrix{ComplexF64}(zeros(N, N))
    ρ0[1, 1] = 1.0

    T = 50.0:50.0:500.0

    β = 1 / (300 * 3.16683e-6) # T = 300K

    Jw = Vector{SpectralDensities.DrudeLorentz}()
    sys_ops = Vector{Matrix{ComplexF64}}()
    for (j, (λ, γ)) in enumerate(zip(λs, γs))
        push!(Jw, SpectralDensities.DrudeLorentz(; λ, γ, Δs=1.0))
        op = zeros(N, N)
        op[j, j] = 1.0
        push!(sys_ops, op)
    end

    ts, ρs = HEOM.propagate(;
                         Hamiltonian=H0,
                         Jw=Jw,
                         β,
                         ρ0,
                         dt,
                         ntimes=nsteps,
                         sys_ops,
                         Lmax=Lmax,
                         num_modes=num_modes
                        )
    
    open("HEOM-SF-populations-Troisi-Rubrene.txt", "w") do io
        pops = [real.(ρs[:, i, i]) for i in 1:N]
        tpops = [ts pops...]
        writedlm(io, tpops, ' ')
    end


end

SingletFissionRedfield()

SingletFissionHEOM()

