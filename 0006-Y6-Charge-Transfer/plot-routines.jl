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


function UpconversionPlotDimer(fname)
    # Plot populations

    dlm = readdlm(fname, Float64)
    
    #boltzdlm = readdlm("../boltzmann-populations.stdout", Float64)
    
    #bXT1 = boltzdlm[1, 1]
    #bXT2 = boltzdlm[2, 1]
    #bCT1 = boltzdlm[3, 1]
    #bCT2 = boltzdlm[4, 1]
    #bTT1 = boltzdlm[5, 1]
    

    #tags = split(fname, "-")
    #s = parse(Float64, tags[6])
    #V = parse(Float64, tags[8])

    t = dlm[:, 1].*au2fs
    
    ρs = []

    N = size(dlm)[2] - 1

    nsteps = size(dlm)[1]-1
    ρXT1 = dlm[:, 2]
    ρXT2 = dlm[:, 3]
    ρCT1 = dlm[:, 4]
    ρCT2 = dlm[:, 5]
    ρTT1 = dlm[:, 6]
    #ρTT2 = dlm[:, 7]
    #plot!(t[2:nsteps], ρTT[2:nsteps], label="TT")
    #plot!(t[2:nsteps], ρXT1[2:nsteps]+ρXT2[2:nsteps], label="XT")
    #plot!(t[2:nsteps], ρCT1[2:nsteps]+ρCT2[2:nsteps], label="CT")
    
    @gp "set key left"
    @gp :- "set title 'Y6 Dimer HEOM" #(Free Parameters - Vct = $V meV, SoC = $s meV)'"
    @gp :- "set xlabel 't(fs)'"
    @gp :- "set ylabel 'Population'"
   
    @gp :- "set yrange[-0.1:1.0]"
    @gp :- "set xrange[0:200]"

    @gp :- t[2:nsteps] ρTT1[2:nsteps] "w l tit 'TT1' dt 1 lw 2 lc rgb 'red' "
   # @gp :- t[2:nsteps] ρTT2[2:nsteps] "w l tit 'TT2' dt 1 lw 2 lc rgb 'red' "
    @gp :- t[2:nsteps] ρXT1[2:nsteps] "w l tit 'XT1' dt 1 lw 2 lc rgb 'cyan' "
    @gp :- t[2:nsteps] ρXT2[2:nsteps] "w l tit 'XT2' dt 1 lw 2 lc rgb 'blue' "
    @gp :- t[2:nsteps] ρCT1[2:nsteps] "w l tit 'CT1' dt 1 lw 2 lc rgb '#74C476' "
    @gp :- t[2:nsteps] ρCT2[2:nsteps] "w l tit 'CT2' dt 1 lw 2 lc rgb '#238B45' "
   
#=
    @gp :- "plot $bTT1 w p notitle lc rgb 'red' pointtype 0"
    @gp :- "plot $bCT1 w p notitle lc rgb '#74C476' pointtype 0"
    @gp :- "plot $bCT2 w p notitle lc rgb '#238B45' pointtype 0"
    @gp :- "plot $bXT1 w p notitle lc rgb 'cyan' pointtype 0"
    @gp :- "plot $bXT2 w p notitle lc rgb 'blue' pointtype 0"
=#
    Gnuplot.save("Results/Y6-HEOM-Dimer-Populations-zoom.png", term="pngcairo size 550,350 fontscale 0.8")

    #savefig("Dimer-TTM-Populations.png")
end


function UpconversionPlotTrimer(fname)
    # Plot populations

    dlm = readdlm(fname, Float64)
    
    t = dlm[:, 1].*au2fs
    
    ρs = []

    N = size(dlm)[2] - 1

    nsteps = size(dlm)[1]-1
    ρXT1 = dlm[:, 2]
    ρXT2 = dlm[:, 3]
    ρXT3 = dlm[:, 4]
    ρCT1 = dlm[:, 5]
    ρCT2 = dlm[:, 6]
    ρCT3 = dlm[:, 7]
    ρCT4 = dlm[:, 8]
    ρCT5 = dlm[:, 9]
    ρCT6 = dlm[:, 10]

    ρXT = ρXT1 + ρXT2 + ρXT3
    ρCT = ρCT1 + ρCT2 + ρCT3 + ρCT4 + ρCT5 + ρCT6

    #plot!(t[2:nsteps], ρTT[2:nsteps], label="TT")
    #plot!(t[2:nsteps], ρXT1[2:nsteps]+ρXT2[2:nsteps], label="XT")
    #plot!(t[2:nsteps], ρCT1[2:nsteps]+ρCT2[2:nsteps], label="CT")
    
    @gp "set key left"
    @gp :- "set title 'Y6 Trimer Redfield'"
    @gp :- "set xlabel 't(fs)'"
    @gp :- "set ylabel 'Population'"
    
    #@gp :- t[2:nsteps] ρTT[2:nsteps] "w l tit 'TT' dt 1 lw 2 lc rgb 'red' "
    @gp :- t[2:nsteps] ρXT1[2:nsteps] "w l tit 'XT1' dt 1 lw 2 lc rgb 'cyan' "
    @gp :- t[2:nsteps] ρXT2[2:nsteps] "w l tit 'XT2' dt 1 lw 2 lc rgb 'blue' "
    @gp :- t[2:nsteps] ρXT3[2:nsteps] "w l tit 'XT3' dt 1 lw 2 lc rgb 'blue' "
    @gp :- t[2:nsteps] ρCT1[2:nsteps] "w l tit 'CT1' dt 1 lw 2 lc rgb '#74C476' "
    @gp :- t[2:nsteps] ρCT2[2:nsteps] "w l tit 'CT2' dt 1 lw 2 lc rgb '#238B45' "
    @gp :- t[2:nsteps] ρCT3[2:nsteps] "w l tit 'CT3' dt 1 lw 2 lc rgb '#238B45' "
    @gp :- t[2:nsteps] ρCT4[2:nsteps] "w l tit 'CT4' dt 1 lw 2 lc rgb '#238B45' "
    @gp :- t[2:nsteps] ρCT5[2:nsteps] "w l tit 'CT5' dt 1 lw 2 lc rgb '#238B45' "
    @gp :- t[2:nsteps] ρCT6[2:nsteps] "w l tit 'CT6' dt 1 lw 2 lc rgb '#238B45' "
    
    Gnuplot.save("Results/all-populations-trimer-brme.png", term="pngcairo size 550,350 fontscale 0.8")

    @gp "set key left"
    @gp :- "set title 'Y6 Trimer Redfield'"
    @gp :- "set xlabel 't(fs)'"
    @gp :- "set ylabel 'Population'"
    
    #@gp :- t[2:nsteps] ρTT[2:nsteps] "w l tit 'TT' dt 1 lw 2 lc rgb 'red' "
    @gp :- t[2:nsteps] ρXT[2:nsteps] "w l tit 'XT' dt 1 lw 2 lc rgb 'cyan' "
    @gp :- t[2:nsteps] ρCT[2:nsteps] "w l tit 'CT' dt 1 lw 2 lc rgb 'blue' "
    
    Gnuplot.save("Results/cum-population-trimer-brme.png", term="pngcairo size 550,350 fontscale 0.8")
    #savefig("Dimer-TTM-Populations.png")
end

function UpconversionYieldDimer()
    files = readdir("stdout-files5")    
    s = zeros(500)
    v = zeros(500)
    yield = zeros(500)

    S = zeros(500)
    V = zeros(500)

    for (i, file) in enumerate(files)
        tags = split(file, "-")
        s[i] = parse(Float64, tags[5]) # SoC in invcm
        v[i] = parse(Float64, tags[7]) # VCT-TT in meV
        
        fpath = "stdout-files5/" * file
        
        dlm = readdlm(fpath, Float64)
        yield[i] =  dlm[size(dlm)[1], size(dlm)[2]]
    end
   #=
    @gsp :- "set view map"
    @gsp :- "set dgrid3d"
    @gsp :- "set key off"
    @gsp :- "set title 'Y6 Dimer Triplet Yield Redfield'"
    @gsp :- "set xlabel 'Spin-Orbit Coupling (invcm)'"
    @gsp :- "set ylabel 'CT-TT Coupling (meV)'"
    
    @gsp :-  "set pm3d interpolate 6, 6" 
    # @gp :- s yield "w l tit 'yield' dt 1 lw 2 lc rgb 'cyan' "

    @gsp :- s v yield "w pm3d tit 'yield'"
    
    Gnuplot.save("Results/yield-vs-soc-and-Vct-Proper.png", term="pngcairo size 550,350 fontscale 0.8")
=#
    @gp "set key left"
    @gp :- "set title 'Y6 Redfield Triplet Yield'"
    @gp :- "set xlabel 'Vct (meV)'"
    @gp :- "set ylabel 'Yield'"
    
    @gp :-  v yield "w p tit 'yield' dt 1 lw 2 lc rgb 'cyan' "
    
    Gnuplot.save("Results/yield-vs-Vct.png", term="pngcairo size 550,350 fontscale 0.8")

    #savefig("Dimer-TTM-Populations.png")

    
end

#for fname in readdir("stdout-files5")
#    UpconversionPlotDimer("stdout-files5/" * fname)
#end


UpconversionPlotDimer("upconversion-dimer-populations-heom.stdout")

#UpconversionYieldDimer()
# Used 0.5nm arbitrarily here, need to figure out site distances for real materials
#plot1DHolstein("Results/Ordejon-Rubrene/ordejon1D-populations.stdout", 0.5)
#UpconversionPlotDimer("Results/upconversion-populations-ttm.stdout")
#UpconversionPlotTrimer("Results/upconversion-populations-brme.stdout")


