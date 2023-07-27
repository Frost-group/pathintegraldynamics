## Testing julia codes for nonadiabatic path integral methods


Some interesting things to try :


- [x] Reimplement a spin-boson model in QuantumDynamics

- [ ] Reimplement a Tully model in QuantumDynamics (no clue how; there's no surface hopping method...Hamiltonian not straightforward in this frame).

- [x] Compare QuAPI with SMatPI + QuAPI from PathSum (Matches most when kmax = 1 in QuantumDynamics ; Is PathSum SMatPI + QuAPI non-Markovian??)

- [ ] Compare QuantumDynamics BlipSum and SMatPI + BlipSum from PathSum (Blips slow)

- [ ] Compare QuAPI + TTM and BlipSum + TTM with QuAPI + SMatPI and BlipSum+ SMatPI for same system. (TTM too slow)

- [ ] Photosynthesis systems (LH2 or BChI or something) model hamiltonian + run method.

- [x] FMO complex from Bose JCP paper

- [ ] FMO complex with QuAPI (couldn't get it to work, some debugging issue)

- [ ] FMO in NQCDynamics
