using LinearAlgebra
using OrdinaryDiffEq
import HierarchicalEOM: Boson_DrudeLorentz_Matsubara, Boson_Underdamped_Matsubara, M_Boson, HEOMsolve, BosonBath, Qobj, ket2dm, basis, getRho, getADO, addTerminator

# As from Lucy; should be cross-checked and linked back to NIST
const hbar = 6.62607015e-34 / (2 * pi)
const e = 1.60e-19
const c = 3e8
const invcm_to_eV = 100 * (hbar * 2 * pi * c) / e
const eV_to_natural = e / hbar
const kT = eV_to_natural * 1.38e-23 * 298 / e
const invcm_to_natural = invcm_to_eV * eV_to_natural

const eV_to_invcm = (100 * (hbar * 2 * pi * c) / e)^-1

function GetHam(dimer_num)

    ESFE = eV_to_invcm * [1.947, 1.9065, 1.868, 1.932, 1.976, 1.9135] #Effective FE singlet state energies in eV D1-D6 Y6
    ETFE = eV_to_invcm * [1.48075, 1.478, 1.4655, 1.46375, 1.50725, 1.469] #Effective FE triplet state energies in eV D1-D6 Y6
    ESCT = eV_to_invcm * [1.675, 1.71, 1.656, 1.7425, 1.8145, 1.7785] #Effective CT singlet energies in eV D1-D6 Y6
    ETCT = eV_to_invcm * [1.7365, 1.749, 1.6705, 1.8035, 1.85, 1.799]

    #Effective SOCME in Y6 Dimer between different states in cm^-1 for Y6 Dimers D1-D6 in order
    VSCTTCT = [0.22838165, 0.11302104, 0.08913591, 0.25666594, 0.08708204, 0.05472136] #SOCME CT(1)-CT(3)
    VSCTTFE = [0.88935921, 0.34432261, 0.31354126, 0.49238916, 0.47722765, 0.35850297] #SOCME CT(1)-FE(3)
    VSFETCT = [0.59736616, 0.17427808, 0.2683069, 0.31014608, 0.11665094, 0.19562612] #SOCME FE(1)-CT(3)
    VSFETFE = [0.09, 0.04582576, 0.01, 0.07211103, 0.07348469, 0.14177447] #SOCME FE(1)-FE(3)
    VSGSTFE = [3.14478736, 2.58708263, 2.33851902, 3.28033231, 2.58540237, 2.00456674] #SOCME GS-FE(3)
    VSGSTCT = [0.80187281, 0, 0.12, 0.30066593, 0.0728011, 0.19672316] #SOCME GS-CT(3)

    #Electron and hole hopping matrix elements of Y6 Dimers D1-D6 in order in eV calculated using counterpoise method
    Vel = eV_to_invcm * [0.05642, -0.05579, -0.02782, 0.03268, -0.04601, -0.03092]
    Vh = eV_to_invcm * [-0.05107, 0.02935, 0.03108, 0.0069, -0.01783, -0.01711]

    VCTFE = 2 * (abs(Vel[dimer_num]) + abs(Vh[dimer_num]))

    #Defining the ground state energy of the Y6 Dimer
    E_S0 = 0

    Ham = [
        ESCT[dimer_num] VCTFE VSCTTCT[dimer_num] VSCTTFE[dimer_num] 0
        VCTFE ESFE[dimer_num] VSFETCT[dimer_num] VSFETFE[dimer_num] 0
        VSCTTCT[dimer_num] VSFETCT[dimer_num] ETCT[dimer_num] VCTFE VSGSTCT[dimer_num]
        VSCTTFE[dimer_num] VSFETFE[dimer_num] VCTFE ETFE[dimer_num] VSGSTFE[dimer_num]
        0 0 VSGSTCT[dimer_num] VSGSTFE[dimer_num] E_S0
    ]

    return Ham
end



function RunHEOM(dimer, tier; t=1000e-15)
    Nkeep = 5 #How many states to keep from the Hamiltonian 
    InputHam = invcm_to_natural * GetHam(dimer)
    InputHam = InputHam[1:Nkeep, 1:Nkeep]

    Hsys = Qobj(InputHam)

    #Specify spectral density function
    #I have based this on the partiton between low and high modes from Samuele's paper
    #Didn't like having the high frequency mode at 0.18 eV though so I have moved it while 
    #conserving the HR factor of the mode. 
    #Conserved the total reorganisationenergy by shifting some to the lower frequency mode

    #Specify spectral density function
    lambda_p = eV_to_natural * 0.034#Samuele
    gamma_p = eV_to_natural * 0.05#Arb. I have used approx 2kT, meant to be for the 'classical' modes
    lambda_u = eV_to_natural * 0.052#Samuele
    gamma_u = eV_to_natural * 0.015 #Arb. (can't be too narrow or won't converge)
    peak_u = eV_to_natural * 0.16  #Samuele 
    baths = BosonBath[]
    operators = []

    include_Ubath = 1

    for i = 1:Nkeep
        Q = ket2dm(basis(Nkeep, i - 1))
        push!(operators, Q)
        DL_bath = Boson_DrudeLorentz_Matsubara(Q, lambda_p, gamma_p, kT, 3)
        push!(baths, DL_bath)
        if include_Ubath == 1
            U_bath = Boson_Underdamped_Matsubara(Q, sqrt(2 * lambda_u * peak_u^2),
                gamma_u, peak_u, kT, 1)
            push!(baths, U_bath)
        end
    end

    #Get Eigenstates
    #Each column in EigenSol.vectors is an eigenstate
    EigenSol = eigen(Hsys)
    state_operators = []
    for i = 1:Nkeep
        state = Qobj(EigenSol.vectors[:, i])
        push!(state_operators, ket2dm(state))
    end

    #Start in singlet FE state
    rho0 = operators[2]

    global L = M_Boson(Hsys, tier, baths)
    for i = 1:Nkeep
        if include_Ubath == 0
            global L = addTerminator(L, baths[i])
        else
            global L = addTerminator(L, baths[2*(i-1)+1])
        end
    end

    #Time list; in seconds
    # Nb, if adaptive, these are just the evaluation points of the ODE
    tlist = 0:t/10000:t

    # WORK DONE HERE
    #  dtmax set explicitly to stop solver adapative step being TOO big and causing Int overflow
    sol = HEOMsolve(L, rho0, tlist, dtmax = 1e-13, e_ops=state_operators, alg=ROCK4())

    return sol 
end

if abspath(PROGRAM_FILE) == @__FILE__
    @time sol=RunHEOM(2, 2, t=1e-12)
    
    result = vcat(tlist', real(sol.expect))
    filename = "Y6_D$(dimer)_EigenSolution_0.txt"
    open(filename, "w") do io
        for i in 1:size(result, 1)
            println(io, result[i, :])
        end
    end

end

