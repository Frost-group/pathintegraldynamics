using OrdinaryDiffEq
using LinearAlgebra

const thz2au = 0.0001519828500716
const invcm2au = 4.55633e-6
const au2fs = 0.02418884254
const mev2invcm = 8.066
const mev2au = mev2invcm * invcm2au
const mev2J = 1.6*10^(-19)*10^(-3)
const nm2au = 18.897

#=

#Old code with incorrect rate matrix

function get_Kmatrix(N, H0, λ, β)
    K = zeros(N, N)
    for a = 1:N
        for b = 1:N
            if a != b
                E = H0[a, a] - H0[b, b] + λ     # ΔG + λ term
                K[a, b] = (H0[a, b])^2 * sqrt(π*β/λ) * exp(-β*(E^2/(4*λ))) 
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
    return K
end


struct Params
    H::AbstractMatrix{ComplexF64}
    K
end

function func_Pauli(p, params, t)
    K = deepcopy(params.K)
    dp = K*p
    return dp
end
=#

struct Params
    kf1f2
    kf2f1
    kfc1
    kcf1
    kfc2
    kcf2
    kct
    ktc
    kft1
    ktf1
    kft2
    ktf2
end

function Params(H, β, λ)
    kf1f2 = Marcus(real.(H[1,2]), real.(H[2,2] - H[1,1]), β, λ[1])
    kf2f1 = Marcus(real.(H[1,2]), real.(H[1,1] - H[2,2]), β, λ[1])
    kfc1 = Marcus(real.(H[1,3]), real.(H[3,3] - H[1,1]), β, λ[3])
    kcf1 = Marcus(real.(H[1,3]), real.(H[1,1] - H[3,3]), β, λ[3])
    kfc2 = Marcus(real.(H[2,3]), real.(H[3,3] - H[2,2]), β, λ[3])
    kcf2 = Marcus(real.(H[2,3]), real.(H[2,2] - H[3,3]), β, λ[3])
    kct = Marcus(real.(H[4,5]), real.(H[5,5] - H[4,4]), β, λ[3])
    ktc = Marcus(real.(H[4,5]), real.(H[4,4] - H[5,5]), β, λ[3])
    kft1 = Marcus(real.(H[1,5]), real.(H[5,5] - H[1,1]), β, λ[5])
    ktf1 = Marcus(real.(H[1,5]), real.(H[1,1] - H[5,5]), β, λ[5])
    kft2 = Marcus(real.(H[2,5]), real.(H[5,5] - H[2,2]), β, λ[5])
    ktf2 = Marcus(real.(H[2,5]), real.(H[2,2] - H[5,5]), β, λ[5])

    Params(kf1f2, kf2f1, kfc1, kcf1, kfc2, kcf2, kct, ktc, kft1, ktf1, kft2, ktf2)
end


Marcus(V, ΔE, β, λ) =2*pi*((V^2 * sqrt(β))/sqrt(4*π*λ))*exp(-1*β*(λ+ΔE)^2/(4*λ))
#/(4.136*10^-12*mev2au)



function func_semi2!(du, u, p, t)
    du[1]=-(p.kf1f2 + 2*p.kfc1 + p.kft1)*u[1]+ p.kf2f1*u[2]+p.kcf1*u[3]+p.kcf1*u[4]+p.ktf1*u[5]
    du[2]=p.kf1f2*u[1] -(p.kf2f1 + 2*p.kfc2 + p.kft2)*u[2]+ p.kcf2*u[3]+ p.kcf2*u[4]+ p.ktf2*u[5]
    du[3]=p.kfc1*u[1]+p.kfc2*u[2] -(p.kcf1+p.kcf2+p.kct)*u[3]+ p.ktc*u[5]
    du[4]=p.kfc1*u[1]+p.kfc2*u[2]  -(p.kcf1+p.kcf2 + p.kct)*u[4]+ p.ktc*u[5]
    du[5]=p.kft1*u[1]+p.kft2*u[2]+p.kct*u[3]+p.kct*u[4] -(p.ktf1+p.ktf2 + 2*p.ktc)*u[5]
end


function propagate(N, H0, λ, β, ntimes, dt)
    #K = get_Kmatrix(N, H0, λ, β)
    #println(K.*invtime2au)
    par = Params(H0, β, λ)
    p0 = zeros(N)
    p0[1] = 1.0
    tspan = (0.0, ntimes*dt)
    prob = ODEProblem(func_semi2!, p0, tspan, par, saveat=0.1)
    sol = solve(prob, Tsit5(), saveat=dt)
    ps = zeros(ComplexF64, length(sol.t), N)
    for j = 1:length(sol.t)
        @inbounds ps[j, :] .= sol.u[j]
    end

    return sol.t, ps
end

