using QuantumDynamics
using LinearAlgebra
using DelimitedFiles
using OrdinaryDiffEq
using HierarchicalEOM
using Gnuplot
using ProgressMeter

# --- Constants ---
#   just ported from PRanay's code, need checking
const thz2au = 0.0001519828500716
const invcm2au = 4.55633e-6
const au2fs = 0.02418884254
const mev2invcm = 8.066
const mev2au = mev2invcm * invcm2au
const nm2au = 18.897

# --- Helper Functions ---

"""
    Marcus(V, ΔE, β, λ)

Calculate the Marcus rate constant. Pranay's code... looks a little odd? check!
"""
Marcus(V, ΔE, β, λ) = 2 * pi * ((V^2 * sqrt(β)) / sqrt(4 * π * λ)) * exp(-1 * β * (λ + ΔE)^2 / (4 * λ))

# --- Hamiltonian Models ---
"""
    Y6Dimer

Struct representing the Y6 Dimer system.

Fields:
- `r::Float64`: Distance parameter.
- `H0::Matrix{ComplexF64}`: Hamiltonian matrix.
- `reorg::Vector{Float64}`: Reorganization energies for each site.
- `cutoff::Vector{Float64}`: Cutoff frequencies for the bath.
"""
struct Y6Dimer
    r::Float64
    H0::Matrix{ComplexF64}
    reorg::Vector{Float64}
    cutoff::Vector{Float64}
end

"""
Samuele's paper fit for CT state energy
"""
Ect(r) = (2.19 - 4.959 / (r)) * 1000 # Result in meV, input r in Angstrom

function Y6Dimer(r::Float64, V::Float64, De::Float64, Dh::Float64)
    Ec = Ect(r)

    # Zhenghan singlet values
    Efe1 = 1872.0
    Efe2 = 1886.0 # never used ?!

    # Zhenghan triplet values
    Et1 = 1350.0
    Et2 = 1393.0

    # Best-guess triplet couplings
    Vt = -76.0 # Triplet-Triplet coupling, assumed to be same as singlet-singlet for now

    Vcte = 15.0 # CT-Triplet coupling (same site)
    Vcth = 15.0 # CT-Triplet coupling (different site)
    socs = 10.0  # Singlet-triplet coupling (same site)
    socn = 10.0  # Singlet-triplet coupling (different site)

    H0 = Matrix{ComplexF64}([
        Efe1 V Dh De socs socn
        V Efe1 De Dh socn socs
        Dh De Ec 0.0 Vcth Vcte
        De Dh 0.0 Ec Vcte Vcth
        socs socn Vcth Vcte Et1 Vt
        socn socs Vcte Vcth Vt Et2
    ]) .* mev2au

    # Reorganization energies and cutoffs (converted to atomic units)
    reorg = [157.0, 157.0, 240.0, 240.0, 157.0, 157.0] .* mev2au
    cutoff = repeat([1600 * invcm2au], 6)

    Y6Dimer(r, H0, reorg, cutoff)
end

function Base.show(io::IO, obj::Y6Dimer)
    println(io, "Y6Dimer(r=$(obj.r))")
    println(io, "H0 dims: $(size(obj.H0))")
end

# --- Rate Matrix Solver ---
struct RateParams
    K::Matrix{Float64}
end

function RateParams(dimer::Y6Dimer, β)
    H = dimer.H0
    λ = dimer.reorg
    N = size(H)[1]
    K = zeros(N, N)
    for i in 1:N
        for j in 1:N
            if i != j
                # Use the larger reorganization energy of the two sites involved
                λ_val = max(λ[i], λ[j])

                K[i, j] = Marcus(real(H[i, j]), real(H[j, j] - H[i, i]), β, λ_val)
            end
        end
    end

    for a in 1:N
        s = 0.0
        for b in 1:N
            s += K[b, a]
        end
        K[a, a] = -s
    end

    return RateParams(K)
end

function func_rate!(du, u, p, t)
    du .= p.K * u
end

# --- Dynamics Methods Types ---

abstract type DynamicsMethod end

struct HEOM <: DynamicsMethod
    L::Int
    K::Int
end
HEOM(; L=4, K=2) = HEOM(L, K)

struct Redfield <: DynamicsMethod end

struct TTM <: DynamicsMethod
    rmax::Int
end
TTM(; rmax=15) = TTM(rmax)

struct Rate <: DynamicsMethod end


# --- Universal Interface ---

"""
    setup_simulation(dimer::Y6Dimer)

Helper to extract common simulation parameters.
"""
function setup_simulation(dimer::Y6Dimer)
    H0 = dimer.H0
    λs = dimer.reorg
    γs = dimer.cutoff
    N = size(H0, 1)

    # T = 300K
    β = 1 / (300 * 3.16683e-6)

    ρ0 = Matrix{ComplexF64}(zeros(N, N))
    ρ0[1, 1] = 1.0

    return H0, λs, γs, N, β, ρ0
end

"""
    run_simulation(method::DynamicsMethod, dimer::Y6Dimer; dt, nsteps)

Run dynamics simulation using the specified method type.
"""
function run_simulation(method::HEOM, dimer::Y6Dimer; dt=0.25 / au2fs, nsteps=4000)
    H0, λs, γs, N, β, ρ0 = setup_simulation(dimer)
    L, K = method.L, method.K
    println("Running HEOM (L=$L, K=$K)...")

    JwH = Vector{SpectralDensities.SpectralDensity}()
    sys_ops = Vector{Matrix{ComplexF64}}()
    for (j, (λ, γ)) in enumerate(zip(λs, γs))
        push!(JwH, SpectralDensities.DrudeLorentz(; λ, γ, Δs=1.0))
        op = zeros(N, N)
        op[j, j] = 1.0
        push!(sys_ops, op)
    end

    times, ρs = QuantumDynamics.HEOM.propagate(;
        Hamiltonian=H0,
        ρ0,
        β,
        dt,
        ntimes=nsteps,
        Jw=JwH,
        sys_ops=sys_ops,
        num_modes=K,
        Lmax=L)


    pops = hcat([real.(ρs[:, k, k]) for k in 1:N]...)
    return times, pops
end

function run_simulation(method::Redfield, dimer::Y6Dimer; dt=0.25 / au2fs, nsteps=4000)
    H0, λs, γs, N, β, ρ0 = setup_simulation(dimer)
    println("Running Redfield...")
    λ_av = sum(λs) / length(λs)
    γ_av = sum(γs) / length(γs)
    η = max(2 * λ_av / (β * γ_av^2), 2 * λ_av / (π * γ_av))
    if η > 1
        @info "Warning: η = $η > 1 ; Redfield might be inaccurate here"
    end

    JwH = Vector{SpectralDensities.SpectralDensity}()
    sys_ops = Vector{Matrix{ComplexF64}}()
    for (j, (λ, γ)) in enumerate(zip(λs, γs))
        push!(JwH, SpectralDensities.DrudeLorentz(; λ, γ, Δs=1.0))
        op = zeros(N, N)
        op[j, j] = 1.0
        push!(sys_ops, op)
    end

    times, ρs = BlochRedfield.propagate(;
        Hamiltonian=H0,
        ρ0,
        β,
        dt,
        ntimes=nsteps,
        Jw=JwH,
        sys_ops=sys_ops
    )

    pops = hcat([real.(ρs[:, k, k]) for k in 1:N]...)
    return times, pops
end

function run_simulation(method::TTM, dimer::Y6Dimer; dt=0.25 / au2fs, nsteps=4000)
    H0, λs, γs, N, β, ρ0 = setup_simulation(dimer)
    rmax = method.rmax
    println("Running TTM (rmax=$rmax)...")

    Jw = [SpectralDensities.DrudeLorentz(; λ=λs[i], γ=γs[i], Δs=1.0) for i in 1:N]

    fbU = Propagators.calculate_bare_propagators(; Hamiltonian=H0, dt=dt, ntimes=nsteps + 1)

    svec = Matrix{Float64}(I, N, N)

    times, ρs = QuantumDynamics.TTM.propagate(;
        fbU=fbU,
        Jw=Jw,
        β=β,
        ρ0=ρ0,
        dt=dt,
        ntimes=nsteps,
        rmax=rmax,
        svec=svec,
        extraargs=TEMPO.TEMPOArgs(),
        path_integral_routine=TEMPO.build_augmented_propagator,
        verbose=true
    )

    pops = hcat([real.(ρs[:, k, k]) for k in 1:N]...)
    return times, pops
end

function run_simulation(method::Rate, dimer::Y6Dimer; dt=0.25 / au2fs, nsteps=4000)
    H0, λs, γs, N, β, ρ0 = setup_simulation(dimer)
    println("Running Rate Matrix...")
    par = RateParams(dimer, β)
    p0 = zeros(Float64, N)
    p0[1] = 1.0
    tspan = (0.0, nsteps * dt)

    prob = ODEProblem(func_rate!, p0, tspan, par)
    sol = solve(prob, Tsit5(), saveat=dt, progress=true)

    times = sol.t
    pops = zeros(length(times), N)
    for t_idx in 1:length(times)
        pops[t_idx, :] = sol.u[t_idx]
    end

    return times, pops
end

# --- I/O and Plotting Interface ---

"""
    save_results(filename, times, populations)

Save simulation results to a delimited file.
"""
function save_results(filename, times, populations)
    open(filename, "w") do io
        writedlm(io, [times populations], ' ')
    end
    println("Saved results to $filename")
end

"""
    load_results(filename)

Load simulation results from a delimited file.
Returns (times, populations).
"""
function load_results(filename)
    data = readdlm(filename, Float64)
    times = data[:, 1]
    populations = data[:, 2:end]
    return times, populations
end

"""
    plot_results(times, populations; title="Population Dynamics", filename=nothing, show=true)

Plot populations using Gnuplot.
If `filename` is provided, saves the plot to that file (e.g. .png).
"""
function plot_results(times, populations; title="Population Dynamics", filename=nothing, show=true)
    # Assumes populations is Size (Time x N)

    times_fs = times .* au2fs

    @gp "set key right"
    @gp :- "set title '$title'"
    @gp :- "set xlabel 't (fs)'"
    @gp :- "set ylabel 'Population'"

    # Define colors/names based on Y6 system knowledge (6 states)
    # 1,2: XT; 3,4: CT; 5,6: TT
    labels = ["XT1", "XT2", "CT1", "CT2", "TT1", "TT2"]
    colors = ["cyan", "blue", "#74C476", "#238B45", "red", "pink"]

    N = size(populations, 2)

    for i in 1:N
        lbl = (i <= length(labels)) ? labels[i] : "State $i"
        col = (i <= length(colors)) ? colors[i] : "black"
        @gp :- times_fs populations[:, i] "w l tit '$lbl' dt 1 lw 2 lc rgb '$col'"
    end

    if !isnothing(filename)
        Gnuplot.save(term="pngcairo size 800,600 fontscale 1.0", filename)
        println("Saved plot to $filename")
    end
end

"""
    plot_from_file(filename; output_file=nothing)

Load results from `filename` and plot them.
"""
function plot_from_file(filename; output_file=nothing)
    times, pops = load_results(filename)
    if isnothing(output_file)
        output_file = replace(filename, ".txt" => ".png")
    end
    plot_results(times, pops; title="Plot from $filename", filename=output_file)
end

function main()
    # Samuele parameters for all 8 dimers; Pranay extractions (?)
    r_vals = [9.29, 13.62, 13.84, 15.49, 18.15, 18.33, 15.44, 15.97]
    V_vals = [-76.0, 5.1, -6.1, -78.7, 50.8, 56.5, -11.5, -9.0]
    De_vals = [72.0, -55.4, -15.0, 53.1, -68.9, -47.9, 0.0, 0.0]
    Dh_vals = [55.7, -45.5, -27.3, -10.9, 24.3, 33.5, 0.0, 0.0]

    methods_to_run = [Rate(), Redfield(), HEOM(), TTM()]

    # Simulation parameters
    # dt = 0.25 fs
    # time = 200 fs
    # nsteps = time / dt = 800
    time = 1
    dt_fs = 0.25
    nsteps = Int(time / dt_fs)
    dt_au = dt_fs / au2fs

    mkpath("plots")
    # Iterate through target dimers
    for i in [2] # list of Samuele dimers
        println("\n=== Simulating Dimer $i ===")
        dimer = Y6Dimer(r_vals[i], V_vals[i], De_vals[i], Dh_vals[i])
        println("Parameters: $dimer")

        for method in methods_to_run
            try
                println("\n--- $(typeof(method)) ---")
                times, pops = @time run_simulation(method, dimer; dt=dt_au, nsteps=nsteps)

                # Output handling
                method_name = string(typeof(method))
                outfile = joinpath("plots", "d$(i)_$(method_name).txt")
                plotfile = joinpath("plots", "d$(i)_$(method_name).png")

                # Save
                save_results(outfile, times, pops)

                # Plot
                plot_results(times, pops;
                    title="Dimer $i ($method_name)",
                    filename=plotfile)

            catch e
                @error "Failed to run $method for Dimer $i: $(string(e))"
                showerror(stdout, e, catch_backtrace())
            end
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
