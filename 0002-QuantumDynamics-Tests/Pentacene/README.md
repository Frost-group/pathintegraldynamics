# Pentacene Dimer Example (Borelli et al. J. Phys. Chem. B 2021, 125, 5397−5407)

Recreated the pentacene dimer example from the JPC B 2021 paper using TTM (Transfer Tensor Method) and QCPI (Quantum-Classical Path Integral) as implemented in QuantumDynamics.jl.


Notes :

* The plots obtained are significantly more oscillatory than the ones obtained in the original paper. 

* The damping nature appears to be similar but the jagged bits of the curve aren't obtained and it appears to be a very smooth damped oscillator. 

Updated Notes (10 Oct) :

* Better results with QCPI on longer time-step.

* Doesn't converge as well and starts to explode after allowing it to propagate for a while.


