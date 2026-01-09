using OrdinaryDiffEq
using LinearAlgebra

const thz2au = 0.0001519828500716
const invcm2au = 4.55633e-6
const au2fs = 0.02418884254
const mev2invcm = 8.066
const mev2au = mev2invcm * invcm2au
const mev2J = 1.6*10^(-19)*10^(-3)
const nm2au = 18.897

# Doesn't work yet, haven't diagnosed it entirely


Marcus(V, ΔE, β, λ) =2*pi*((V^2 * sqrt(β))/sqrt(4*π*λ))*exp(-1*β*(λ+ΔE)^2/(4*λ))


struct Params
    K
end

function Params(H, β, λ)
    N = size(H)[1]
    K = zeros(N, N)
    for i in 1:N
        for j in 1:N
            l = 3
            if (i == 1 || i == 2) && (j == 1 || j == 2)
                l = 1
            end
            if (i == 5 || i == 6) || (j == 5 || j == 6)
                l = 5
            end
            if (i == 3 || i == 4) || (j == 3 || j == 4)
                l = 3
            end    
            if (i != j)
                K[i, j] = Marcus(real.(H[i, j]), real.(H[j, j] - H[i, i]), β, λ[l])
            end
        end
    end


    for a = 1:N
        s = 0.0
        for b = 1:N
            s += K[b, a]
        end
        K[a, a] = -s
    end
    

    return Params(K)
end


#/(4.136*10^-12*mev2au)


# Something's wrong I can feel it
function func_semi!(du, u, p, t)
    du = p.K * u
end

function propagate(N, H0, λ, β, ntimes, dt)
    par = Params(H0, β, λ)
    println(par.K)
    p0 = zeros(N)
    p0[1] = 1.0
    tspan = (0.0, ntimes*dt)
    prob = ODEProblem(func_semi!, p0, tspan, par, saveat=0.1)
    sol = solve(prob, Tsit5(), saveat=dt)
    ps = zeros(ComplexF64, length(sol.t), N)
    for j = 1:length(sol.t)
        @inbounds ps[j, :] .= sol.u[j]
    end

    return sol.t, ps
end

