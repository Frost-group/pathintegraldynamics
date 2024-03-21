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
    
    nsteps = 5000

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

    nsteps = 2000
    ρXT1 = dlm[:, 2]
    ρXT2 = dlm[:, 3]
    ρCT1 = dlm[:, 4]
    ρCT2 = dlm[:, 5]
    ρTT = dlm[:, 6]
    #plot!(t[2:nsteps], ρTT[2:nsteps], label="TT")
    #plot!(t[2:nsteps], ρXT1[2:nsteps]+ρXT2[2:nsteps], label="XT")
    #plot!(t[2:nsteps], ρCT1[2:nsteps]+ρCT2[2:nsteps], label="CT")
    
    @gp "set key left"
    @gp :- "set title 'Y6 Dimer Redfield (Free Parameters - Vct = 162meV, SoC = 71.0meV)'"
    @gp :- "set xlabel 't(fs)'"
    @gp :- "set ylabel 'Population'"
    
    @gp :- t[2:nsteps] ρTT[2:nsteps] "w l tit 'TT' dt 1 lw 2 lc rgb 'red' "
    @gp :- t[2:nsteps] ρXT1[2:nsteps] "w l tit 'XT1' dt 1 lw 2 lc rgb 'cyan' "
    @gp :- t[2:nsteps] ρXT2[2:nsteps] "w l tit 'XT2' dt 1 lw 2 lc rgb 'blue' "
    @gp :- t[2:nsteps] ρCT1[2:nsteps] "w l tit 'CT1' dt 1 lw 2 lc rgb '#74C476' "
    @gp :- t[2:nsteps] ρCT2[2:nsteps] "w l tit 'CT2' dt 1 lw 2 lc rgb '#238B45' "
    
    Gnuplot.save("Results/Y6-dimer-upconversion-populations-redfield.png", term="pngcairo size 550,350 fontscale 0.8")

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
    V = zeros(500)
    yield = zeros(500)

    for (i, file) in enumerate(files)
        tags = split(file, "-")
        s[i] = parse(Float64, tags[5]) # SoC in invcm
        V[i] = parse(Float64, tags[7]) # VCT-TT in meV
        
        fpath = "stdout-files5/" * file

        dlm = readdlm(fpath, Float64)
        yield[i] =  dlm[size(dlm)[1], size(dlm)[2]]
    end
   
    @gsp :- "set view map"
    @gsp :- "set dgrid3d"
    @gsp :- "set key off"
    @gsp :- "set title 'Y6 Dimer Triplet Yield Redfield'"
    @gsp :- "set xlabel 'Spin-Orbit Coupling (invcm)'"
    @gsp :- "set ylabel 'CT-TT Coupling (meV)'"
    
    @gsp :-  "set pm3d interpolate 6, 6" 
    # @gp :- s yield "w l tit 'yield' dt 1 lw 2 lc rgb 'cyan' "

    @gsp :- s V yield "w pm3d tit 'yield'"
    
    Gnuplot.save("Results/yield-vs-soc-and-Vct-Proper.png", term="pngcairo size 550,350 fontscale 0.8")

#=    @gp "set key left"
    @gp :- "set title 'Y6 Trimer Redfield'"
    @gp :- "set xlabel 'SOC (meV)'"
    @gp :- "set ylabel 'Yield'"
    
    @gp :-  s yield "w p tit 'yield' dt 1 lw 2 lc rgb 'cyan' "
    
    Gnuplot.save("Results/yield-vs-SoC.png", term="pngcairo size 550,350 fontscale 0.8")=#
    #savefig("Dimer-TTM-Populations.png")

    
end

UpconversionPlotDimer("stdout-files5/upconversion-populations-brme-s-71.0-V-162.0-t-2000.0.stdout")

#UpconversionYieldDimer()
# Used 0.5nm arbitrarily here, need to figure out site distances for real materials
#plot1DHolstein("Results/Ordejon-Rubrene/ordejon1D-populations.stdout", 0.5)
#UpconversionPlotDimer("Results/upconversion-populations-ttm.stdout")
#UpconversionPlotTrimer("Results/upconversion-populations-brme.stdout")


