const thz2au = 0.0001519828500716
const invcm2au = 4.55633e-6
const au2fs = 0.02418884254
const mev2invcm = 8.066
const mev2au = mev2invcm * invcm2au
const nm2au = 18.897



function Y6UpconversionDimerHamiltonian(socs::Float64, socn::Float64, Vcte::Float64, Vcth::Float64)
    
    # Params from Samuele paper

    # D1 dimer dynamics - FE vs CT character over time.
    
    reorg = [157.0, 157.0, 240.0, 240.0, 157.0, 157.0]*mev2au

    cutoff = repeat([1600 * invcm2au], 5)

    #Efe = 2046.0
    Ect(r) = (2.19 -  4.959/(r))*1000  # Enter r in angstrom ; best fit equation for Ect
    Dh = 55.7
    De = 72.0
    V = -76.0
    
    Vt = -76.0 # Triplet-Triplet coupling, assumed to be same as singlet-singlet for now
   
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
        Efe1 V Dh De socs socn 
        V Efe1 De Dh socn socs 
        Dh De Ec 0.0 Vcth Vcte 
        De Dh 0.0 Ec Vcte Vcth
        socs socn Vcth Vcte Et1 Vt
        socn socs Vcte Vcth Vt Et2
    ]) * mev2au 


    return reorg, cutoff, H0
end


function Y6UpconversionDimerSink(socs::Float64, socn::Float64, Vcte::Float64, Vcth::Float64)
    
    # Params from Samuele paper

    # D1 dimer dynamics - FE vs CT character over time.
    
    reorg = [157.0, 157.0, 240.0, 240.0, 157.0]*mev2au

    cutoff = repeat([1600 * invcm2au], 5)

    #Efe = 2046.0
    Ect(r) = (2.19 -  4.959/(r))*1000  # Enter r in angstrom ; best fit equation for Ect
    Dh = 55.7
    De = 72.0
    V = -76.0
    
    Vt = -76.0 # Triplet-Triplet coupling, assumed to be same as singlet-singlet for now
   
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

    Egs = 0.0 # Energy of sink state
    
    Vfg = 100.0 # Strongly couples??
    Vcg = 100.0
    Vtg = 100.0

    H0 = Matrix{ComplexF64}([
        Efe1 V Dh De socs socn Vfg
        V Efe1 De Dh socn socs Vfg
        Dh De Ec 0.0 Vcth Vcte Vfg
        De Dh 0.0 Ec Vcte Vcth Vfg
        socs socn Vcth Vcte Et1 Vt Vfg
        socn socs Vcte Vcth Vt Et2 Vfg
        Vfg Vfg Vcg Vcg Vtg Vtg Egs
    ]) * mev2au 


    return reorg, cutoff, H0
end

function Y6UpconversionTrimerHamiltonian()
    
    # D1, D3, D8 trimer dynamics - FE vs CT character over time.

    λ = [repeat([157*mev2au], 3)..., repeat([240*mev2au], 6)...]

    γ = repeat([1600 * invcm2au], 9)

    Ef1 = 2046.0 # eV
    Ef3 = 2046.0
    Ef8 = 2046.0

    Ect(r) = (2.19 -  4.959/(r))*1000  # Enter r in angstrom ; best fit equation for Ect
   
    r1 = 9.29
    r3 = 13.84
    r8 = 15.97
    Dh1 = 55.7
    Dh3 = -27.3
    Dh8 = 0.0
    De1 = 72.0
    De3 = -15.0
    De8 = 0.0
   
    V1 = -76.0
    V3 = -6.1
    V8 = -9.0
   
    Ec1 = Ect(r1) # 1656.2 eV
    Ec3 = Ect(r3) # 1831.69 eV
    Ec8 = Ect(r8) # 1879.48 eV
   
    te1 = 66.1
    te3 = -27.0
    te8 = 0.0
    th1 = 49.6
    th3 = -11.3
    th8 = 0.0

    
    H0 = Matrix{ComplexF64}([
         Ef1 V1 V8 Dh1 Dh8 De1 0.0 De8 0.0
         V1 Ef3 V3 De1 0.0 Dh1 Dh3 0.0 De3
         V8 V3 Ef8 0.0 De8 0.0 De3 Dh8 Dh3
         Dh1 De1 0.0 Ec1 th3 0.0 0.0 0.0 te8
         Dh8 0.0 De8 th3 Ec8 0.0 te1 0.0 0.0  
         De1 Dh1 0.0 0.0 0.0 Ec1 th8 te3 0.0
         0.0 Dh3 De3 0.0 te1 th8 Ec3 0.0 0.0
         De8 0.0 Dh8 0.0 0.0 te3 0.0 Ec8 th1
         0.0 De3 Dh3 te8 0.0 0.0 0.0 th1 Ec3
        ])
    
    return λ, γ, H0

end


