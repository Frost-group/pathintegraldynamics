using Plots
using DelimitedFiles

function conv()
    dlms = []
    rmax = []
    for fname in readdir()
        if fname[1:4] == "CONV"
            push!(dlms, readdlm(fname, Float64))
            push!(rmax, parse(Float64, fname[30:31]))
        end 

    end
    


    for (i, dlm) in enumerate(dlms)
        r = rmax[i]
        N = size(dlm)[2] - 1
        nsteps = size(dlm)[1]-1
        t = dlm[:,1]
        ρ4 = dlm[:, 4]
        plot!(t[2:nsteps], ρ4[2:nsteps], label="site 4 population, rmax=$r",legend=:topleft, xscale=:log10)
    end

    savefig("convergence-plot.png")

end

conv()


