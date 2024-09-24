using Gnuplot
using DelimitedFiles
#include("long-ratematrix.jl")
include("ratematrix-solver.jl")
include("hamiltonian-models.jl")

function UpconvertRateMatrix()
    λs, γ, H0 = Y6UpconversionDimerHamiltonian(71.0, 71.0, 162.0, 162.0)

    N = 6
    
    β = 1052.0 # 300K
    ts, ps = propagate(N, H0, λs, β, 4000, 0.25/au2fs)
    
    ts = ts.*au2fs

    open("upconversion-marcus-populations.stdout", "w") do io
       pops = [real.(ps[:, i]) for i in 1:N]
       tpops = [ts pops...]
       writedlm(io, tpops, ' ')
   end    


    println(size(ts))
    println(size(ps[:, 1]))
    #println(ps)
    #println(ts)
    @gp "set key left"
    @gp :- "set title 'Y6 Dimer Marcus" #(Free Parameters - Vct = $V meV, SoC = $s meV)'"
    @gp :- "set xlabel 't(fs)'"
    @gp :- "set ylabel 'Population'"
    @gp :- "set xrange [0:100]"

    @gp :- ts real.(ps[:, 5]) "w l tit 'TT1' dt 1 lw 2 lc rgb 'red' "
    @gp :- ts real.(ps[:, 6]) "w l tit 'TT2' dt 1 lw 2 lc rgb 'pink' "
    @gp :- ts real.(ps[:, 1]) "w l tit 'XT1' dt 1 lw 2 lc rgb 'cyan' "
    @gp :- ts real.(ps[:, 2]) "w l tit 'XT2' dt 1 lw 2 lc rgb 'blue' "
    @gp :- ts real.(ps[:, 3]) "w l tit 'CT1' dt 1 lw 2 lc rgb '#74C476' "
    @gp :- ts real.(ps[:, 4]) "w l tit 'CT2' dt 1 lw 2 lc rgb '#238B45' "
  
   Gnuplot.save("Y6-dimer-ratematrix.png", term="pngcairo size 550,350 fontscale 0.8")

end



