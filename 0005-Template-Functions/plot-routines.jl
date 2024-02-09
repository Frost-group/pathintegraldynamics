# Simple plotting routines for various systems non-adiabatic dynamics.

using DelimitedFiles
using Plots

const thz2au = 0.0001519828500716
const invcm2au = 4.55633e-6
const au2fs = 0.02418884254
const mev2invcm = 8.066
const mev2au = mev2invcm * invcm2au
const nm2au = 18.897


"""
plot1DHolstein(N, a)

Plot out the populations, MSD and dMSD/dt for an N-site 1D lattice with spacing a (nm).

"""

function plot1DHolstein(N, a)
    # Plot populations

    dlm = readdlm("polaron-1d-pops-rmax10.txt", Float64)
    
    t = dlm[:, 1].*au2fs
    
    ρs = []

    nsteps = size(dlm)[1]-1

    for i in 2:size(dlm)[2]
        ρ = dlm[:, i]
        push!(ρs, ρ)
        plot!(t[2:nsteps], ρ[2:nsteps], label="site $(i-1)", xscale=:log10)
    end
    savefig("Populations.png")
    
    MSD = []
    for i in 1:(nsteps+1)
        s= 0.0
        for j in 1:N
            s += ρs[j][i]*((j-1)^2)*(a * nm2au)^2
        end
        push!(MSD, s)
    end

    plot(t[2:nsteps], MSD[2:nsteps], xscale=:log10)

    savefig("MSD.png")

    dMSD = [(MSD[i] - MSD[i-1])/(t[i] - t[i-1]) for i in 2:(nsteps+1)]

    plot(t[2:nsteps], dMSD[2:nsteps].*(5), xscale=:log10)

    savefig("dMSD_dt.png")



end

plot1DHolstein(10, 0.5)


