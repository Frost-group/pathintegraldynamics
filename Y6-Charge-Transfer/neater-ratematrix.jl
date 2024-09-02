using OrdinaryDiffEq
using LinearAlgebra

const thz2au = 0.0001519828500716
const invcm2au = 4.55633e-6
const au2fs = 0.02418884254
const mev2invcm = 8.066
const mev2au = mev2invcm * invcm2au
const mev2J = 1.6*10^(-19)*10^(-3)
const nm2au = 18.897


Marcus(V, ΔE, β, λ) =2*pi*((V^2 * sqrt(β))/sqrt(4*π*λ))*exp(-1*β*(λ+ΔE)^2/(4*λ))


struct Params
    Kmatrix
end

function Params(H, β, λ)
    N = size(H)[1]
    K = Matrix{ComplexF64}(zeros(N, N))
    for i in 1:N
        for j in 1:N
            if i != j
                k = (i > j) ? i : j
                K[i, j] = Marcus(real.(H[i, j]), real.(H[j, j] - H[i, i]), β, λ[k])
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


function func_semi!(du, u, p, t)
    du = K*u
end


function func_semi2!(du, u, p, t)
    du[1]=-(p.kf1f2 + p.kf1c1 + p.kf1c2+ p.kf1t1+p.kf1t2)*u[1]+p.kf2f1*u[2]+p.kc1f1*u[3]+p.kc2f1*u[4]+p.kt1f1*u[5]+p.kt2f1*u[5]
    du[2]=p.kf1f2*u[1] -(p.kf2f1 + p.kf2c1+ p.kf2c2 + p.kf2t1+p.kf2t2)*u[2]+ p.kc1f2*u[3]+ p.kc2f2*u[4]+ p.kt1f2*u[5]+ p.kt2f2*u[5]
    du[3]=p.kf1c1*u[1]+p.kf2c1*u[2] -(p.kc1f1+p.kc1f2+p.kc1t1+p.kc1t2)*u[3]+ p.kt1c1*u[5]+p.kt2c1*u[5]
    du[4]=p.kf1c2*u[1]+p.kf2c2*u[2] -(p.kc2f1+p.kc2f2+p.kc2t1+p.kc2t2)*u[4]+ p.kt1c2*u[5]+p.kt2c2*u[5]
    du[5]=p.kf1t1*u[1]+p.kf2t1*u[2]+p.kc1t1*u[3]+p.kc2t1*u[4] -(p.kt1f1+p.kt1f2 +p.kt1c1+p.kt1c2+p.kt1t2)*u[5]+p.kt2t1*u[6]
    du[6]=p.kf1t2*u[1]+p.kf2t2*u[2]+p.kc1t2*u[3]+p.kc2t2*u[4] -(p.kt2f1+p.kt2f2 +p.kt2c1+p.kt2c2+p.kt2t1)*u[6]+p.kt1t2*u[5]
end


function propagate(N, H0, λ, β, ntimes, dt)
    #K = get_Kmatrix(N, H0, λ, β)
    #println(K.*invtime2au)
    par = Params(H0, β, λ)
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

