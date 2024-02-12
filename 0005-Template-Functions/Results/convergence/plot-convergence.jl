using Plots
using DelimitedFiles

function conv()
    fname = "CONVERGENCE-populations-rmax_"
    dlm1 = readdlm(fname*"1.txt", Float64)
    dlm2 = readdlm(fname*"2.txt", Float64)
    dlm3 = readdlm(fname*"3.txt", Float64)
    dlm4 = readdlm(fname*"4.txt", Float64)
    dlm10 = readdlm(fname*"10.txt", Float64)
   
    t1 = dlm1[:, 1]
    ρ1 = dlm1[:, 4]

    t2 = dlm2[:, 1]
    ρ2 = dlm2[:, 4]
   
    t3 = dlm3[:, 1]
    ρ3 = dlm3[:, 4]

    t4 = dlm4[:, 1]
    ρ4 = dlm4[:, 4]
    
    t10 = dlm10[:, 1]
    ρ10 = dlm10[:, 4]

    plot!(t1, ρ1, label="rmax1")
    plot!(t2, ρ2, label="rmax2")
    plot!(t3, ρ3, label="rmax3")
    plot!(t4, ρ4, label="rmax4")
    plot!(t10, ρ10, label="rmax10")
    savefig("convergence-plot.png")
end

conv()


