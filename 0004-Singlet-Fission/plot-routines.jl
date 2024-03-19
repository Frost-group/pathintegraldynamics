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


function SFPlotDimer(fname)
    # Plot populations

    dlm = readdlm(fname, Float64)
    
    t = dlm[:, 1].*au2fs
    
    ρs = []

    N = size(dlm)[2] - 1

    nsteps = size(dlm)[1]-1

    ρXT1 = dlm[:, 2]
    ρXT2 = dlm[:, 3]
    ρCT1 = dlm[:, 4]
    ρCT2 = dlm[:, 5]
    ρTT  = dlm[:, 6]
    #plot!(t[2:nsteps], ρTT[2:nsteps], label="TT")
    #plot!(t[2:nsteps], ρXT1[2:nsteps]+ρXT2[2:nsteps], label="XT")
    #plot!(t[2:nsteps], ρCT1[2:nsteps]+ρCT2[2:nsteps], label="CT")
    
    @gp "set key left"
    @gp :- "set title 'Rubrene Dimer Redfield'"
    @gp :- "set xlabel 't(fs)'"
    @gp :- "set ylabel 'Population'"
    
    #@gp :- t[2:nsteps] ρTT[2:nsteps] "w l tit 'TT' dt 1 lw 2 lc rgb 'red' "
    @gp :- t[2:nsteps] ρXT1[2:nsteps] "w l tit 'XT1' dt 1 lw 2 lc rgb 'cyan' "
    @gp :- t[2:nsteps] ρXT2[2:nsteps] "w l tit 'XT2' dt 1 lw 2 lc rgb 'blue' "
    @gp :- t[2:nsteps] ρCT1[2:nsteps] "w l tit 'CT1' dt 1 lw 2 lc rgb '#74C476' "
    @gp :- t[2:nsteps] ρCT2[2:nsteps] "w l tit 'CT2' dt 1 lw 2 lc rgb '#238B45' "
    @gp :- t[2:nsteps]  ρTT[2:nsteps] "w l tit 'TT' dt 1 lw 2 lc rgb 'red' "
    
    Gnuplot.save("Results/populations-rubrene-troisi-redfield.png", term="pngcairo size 550,350 fontscale 0.8")

    #savefig("Dimer-TTM-Populations.png")
end

SFPlotDimer("Results/redfield-SF-populations-Troisi-Rubrene.stdout")


