# Simple plotting routines for various systems non-adiabatic dynamics.

using DelimitedFiles
#using Plots
using Gnuplot

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

function plot1DHolstein(fname, a)
    # Plot populations

    dlm = readdlm(fname, Float64)
    
    t = dlm[:, 1].*au2fs
    
    ρs = []

    N = size(dlm)[2] - 1

    nsteps = size(dlm)[1]-1
    nsteps = 10000
    @gp "set key left" "set logscale x"
    @gp :- "set title 'Populations'"
    for i in 2:size(dlm)[2]
        ρ = dlm[:, i]
        push!(ρs, ρ)
        #plot!(t[2:nsteps], ρ[2:nsteps], label="site $(i-1)", xscale=:log10)
        @gp :- t[2:nsteps] ρ[2:nsteps] "w l tit 'site $(i-1)' dt 1 lw 2"
    end
    #savefig("Populations.png")

    Gnuplot.save("Populations.png", term="pngcairo size 550,350 fontscale 0.8")
    
    MSD = Vector{Float64}(undef, nsteps+1)
    for i in 1:(nsteps+1)
        s= 0.0
        for j in 1:N
            s += ρs[j][i]*((j-1)^2)*(a * nm2au)^2
        end
        MSD[i] = s
    end

    dMSD = [(MSD[i+1] - MSD[i-1])/(t[i+1] - t[i-1]) for i in 2:(nsteps)]

    @gp "set key left" "set logscale x"
    @gp :- "set title 'MSD'"

    @gp :- t[2:nsteps-1] MSD[2:nsteps-1] "w l tit ' ' dt 1 lw 2"

    Gnuplot.save("MSD.png", term="pngcairo size 550,350 fontscale 0.8")
    #plot(t[2:nsteps], MSD[2:nsteps], xscale=:log10)

    #savefig("MSD.png")


    @gp "set key left" "set logscale x"
    @gp :- "set title 'dMSD/dt'"

    @gp :- t[2:nsteps-1] dMSD[2:nsteps-1] "w l tit 'dMSD/dt' dt 1 lw 2"

    Gnuplot.save("dMSD_dt.png", term="pngcairo size 550,350 fontscale 0.8")
    #plot(t[2:nsteps], MSD[2:nsteps], xscale=:log10)
    #plot(t[2:nsteps-1], dMSD[2:nsteps-1].*(5), xscale=:log10)

   # savefig("dMSD_dt.png")
    
   #=
    # Moving average dMSD_dt
    
    n_bin = 45

    mov = []
    for i in 1:(length(enumerate(MSD))-n_bin)
        sum = 0
        for j in 1:n_bin
            sum += MSD[i+j-1]
        end
        push!(mov, sum/(n_bin))
    end
    
    dMSD_agg = [(MSD[i+1] - MSD[i-1])/(t[i+1] - t[i-1]) for i in 2:(nsteps-n_bin-1)]
    
    #plot(t[2:(nsteps-n_bin-1)], dMSD_agg, xscale=:log10)

    savefig("dMSD_agg.png")
    =#
    
end


function UpconversionPlotDimer(fname)
    # Plot populations

    dlm = readdlm(fname, Float64)
    
    t = dlm[:, 1].*au2fs
    
    ρs = []

    N = size(dlm)[2] - 1

    nsteps = size(dlm)[1]-1

    ρTT = dlm[:, 2]
    ρXT = dlm[:, 3] + dlm[:, 4]
    ρCT = dlm[:, 5] + dlm[:, 6]

    #plot!(t[2:nsteps], ρTT[2:nsteps], label="TT")
    #plot!(t[2:nsteps], ρXT1[2:nsteps]+ρXT2[2:nsteps], label="XT")
    #plot!(t[2:nsteps], ρCT1[2:nsteps]+ρCT2[2:nsteps], label="CT")
    
    @gp "set key left" "set logscale x"
    @gp :- "set title 'Dimer Singlet Fission'"
    
    @gp :- t[2:nsteps] ρTT[2:nsteps] "w l tit 'TT' dt 1 lw 2 lc rgb 'red' "
    @gp :- t[2:nsteps] ρXT[2:nsteps] "w l tit 'XT' dt 1 lw 2 lc rgb 'blue' "
    @gp :- t[2:nsteps] ρCT[2:nsteps] "w l tit 'CT' dt 1 lw 2 lc rgb 'green' "
    
    Gnuplot.save("Dimer-TTM-Populations.png", term="pngcairo size 550,350 fontscale 0.8")

    #savefig("Dimer-TTM-Populations.png")
end

# Used 0.5nm arbitrarily here, need to figure out site distances for real materials
plot1DHolstein("Results/ordejon1D-populations-gamma25.stdout", 0.5)
#UpconversionPlotDimer("dimer-populations.stdout")


