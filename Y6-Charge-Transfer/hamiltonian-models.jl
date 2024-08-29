const thz2au = 0.0001519828500716
const invcm2au = 4.55633e-6
const au2fs = 0.02418884254
const mev2invcm = 8.066
const mev2au = mev2invcm * invcm2au
const nm2au = 18.897



function Y6UpconversionDimerHamiltonian(socs::Float64, socn::Float64, Vcte::Float64, Vcth::Float64)
    
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


    H0 = Matrix{ComplexF64}([
        Efe1 V Dh De socs socn Vfg
        V Efe1 De Dh socn socs Vfg
        Dh De Ec 0.0 Vcth Vcte Vfg
        De Dh 0.0 Ec Vcte Vcth Vfg
        socs socn Vcth Vcte Et1 Vt Vfg
        socn socs Vcte Vcth Vt Et2 Vfg
        Vfg Vfg Vfg Vfg Vfg Vfg Egs
    ]) * mev2au 


    return reorg, cutoff, H0
end
