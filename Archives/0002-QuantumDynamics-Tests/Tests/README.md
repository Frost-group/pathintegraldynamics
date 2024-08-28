## Testing julia codes for nonadiabatic path integral methods


Some interesting things to try :


- [x] Reimplement a spin-boson model in QuantumDynamics

- [ ] Reimplement a Tully model in QuantumDynamics (no clue how; there's no surface hopping method...Hamiltonian not straightforward in this frame).

- [x] Compare QuAPI with SMatPI + QuAPI from PathSum (Matches most when kmax = 1 in QuantumDynamics ; Is PathSum SMatPI + QuAPI non-Markovian??)

- [ ] Compare QuantumDynamics BlipSum and SMatPI + BlipSum from PathSum (Blips slow)

- [x] Compare QuAPI + TTM and BlipSum + TTM with QuAPI + SMatPI and BlipSum+ SMatPI for same system.

- [x] FMO complex from Bose JCP paper

- [x] FMO complex with QuAPI (works with TTM but not with plain QuAPI)

- [ ] FMO in NQCDynamics

