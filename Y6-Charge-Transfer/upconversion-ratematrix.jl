using Gnuplot

const thz2au = 0.0001519828500716
const invcm2au = 4.55633e-6
const au2fs = 0.02418884254
const mev2invcm = 8.066
const mev2au = mev2invcm * invcm2au
const nm2au = 18.897


include("ratematrix-solver.jl")

function Y6UpconversionDimerHamiltonian(soc::Float64, Vct::Float64)
    
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

    
    # Dimer Hamiltonian including triplet with placeholder SoC
    
    # Zhenghan singlet values

    Efe1 = 1872.0
    Efe2 = 1886.0

    # Zhenghan triplet values
    Et1 = 1350.0
    Et2 = 1393.0
    
    #Ett = Et1 + Et2

    Ett = 1350.0

    H0 = Matrix{ComplexF64}([
        Efe1 V Dh De soc 
        V Efe1 De Dh soc 
        Dh De Ec 0.0 Vct 
        De Dh 0.0 Ec Vct 
        soc soc Vct Vct Ett
    ]) * mev2au 


    return reorg, cutoff, H0
end



function UpconvertRateMatrix()
    λs, γ, H0 = Y6UpconversionDimerHamiltonian(162.0, 71.0)

    N = 5
    
    β = 1052.0 # 300K
    ts, ps = propagate(N, H0, λs, β, 4000, 0.25/au2fs)
    
    ts = ts.*au2fs
    
    println(size(ts))
    println(size(ps[:, 1]))
    #println(ps)
    #println(ts)
    @gp "set key left"
    @gp :- "set title 'Y6 Dimer Redfield" #(Free Parameters - Vct = $V meV, SoC = $s meV)'"
    @gp :- "set xlabel 't(fs)'"
    @gp :- "set ylabel 'Population'"
    @gp :- "set xrange [0:100]"

    @gp :- ts real.(ps[:, 5]) "w l tit 'TT1' dt 1 lw 2 lc rgb 'red' "
   # @gp :- t[2:nsteps] ρTT2[2:nsteps] "w l tit 'TT2' dt 1 lw 2 lc rgb 'red' "
   @gp :- ts real.(ps[:, 1]) "w l tit 'XT1' dt 1 lw 2 lc rgb 'cyan' "
   @gp :- ts real.(ps[:, 2]) "w l tit 'XT2' dt 1 lw 2 lc rgb 'blue' "
   @gp :- ts real.(ps[:, 3]) "w l tit 'CT1' dt 1 lw 2 lc rgb '#74C476' "
   @gp :- ts real.(ps[:, 4]) "w l tit 'CT2' dt 1 lw 2 lc rgb '#238B45' "
  
   Gnuplot.save("Y6-dimer-ratematrix.png", term="pngcairo size 550,350 fontscale 0.8")

end


UpconvertRateMatrix()
