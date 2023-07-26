## Testing julia codes for nonadiabatic path integral methods


Some interesting things to try :


- [x] Reimplement a spin-boson model

- [ ] Reimplement a Tully model (no clue how; there's no surface hopping method...Hamiltonian not straightforward in this frame).

- [x] Compare QuAPI with SMatPI + QuAPI from PathSum (Matches most when kmax = 1 in QuantumDynamics ; Is PathSum SMatPI + QuAPI non-Markovian??)

- [ ] Compare QuantumDynamics BlipSum and SMatPI + BlipSum from PathSum (Blips slow)

- [ ] Compare QuAPI + TTM and BlipSum + TTM with QuAPI + SMatPI and BlipSum+ SMatPI for same system. (TTM too slow)

- [ ] Photosynthesis PS-2 systems.

- [x] FMO complex from Bose JCP paper

- [ ] FMO complex with QuAPI (couldn't get it to work, some debugging issue)
