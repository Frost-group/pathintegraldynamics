using OrdinaryDiffEq
using LinearAlgebra

function get_Kmatrix(N, H0, λ, β)
    K = zeros(N, N)
    for a = 1:N
        for b = 1:N
            if a != b
                E = H0[a, a] - H0[b, b]
                K[a, b] = (H0[a, b])^2 * sqrt(π*β/λ) * exp(-β*(E + λ)^2/(4*λ))
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


function propagate(N, H0, λ, β, ntimes, dt)
    K = get_Kmatrix(N, H0, λ, β)
    par = Params(H0, K)
    p0 = zeros(N)
    p0[1] = 1.0
    tspan = (0.0, ntimes*dt)
    prob = ODEProblem(func_Pauli, p0, tspan, par)
    sol = solve(prob, Tsit5(), saveat=dt)
    ps = zeros(ComplexF64, length(sol.t), N)
    for j = 1:length(sol.t)
        @inbounds ps[j, :] .= sol.u[j]
    end

    return sol.t, ps
end

