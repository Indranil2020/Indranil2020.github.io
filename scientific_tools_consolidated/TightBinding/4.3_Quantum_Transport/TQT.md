# TQT (Twente Quantum Transport)

## Official Resources
- **Repository**: https://github.com/TwenteQT/TwenteQuantumTransport
- **License**: GNU General Public License v3.0

## Overview
**TQT** ( Twente Quantum Transport) is a versatile, high-performance Fortran code for simulating **spin-dependent electron transport** in nanoelectronics and spintronic devices. It specifically targets realistic material systems with disorder (chemical, thermal, magnetic) by employing a **Scattering Matrix** approach combined with **Green's functions**. TQT supports various Hamiltonian types, including tight-binding (Slater-Koster) and Muffin-Tin Orbitals (LMTO/EMTO), making it suitable for both model studies and first-principles-based transport calculations in large supercells.

**Scientific domain**: Spintronics, Quantum Transport, Disordered Systems
**Target user community**: Researchers in spintronics, nanomagnetism, and device physics

## Theoretical Methods
- **Landauer-Büttiker Formalism**: Calculates conductance and transmission probabilities from the scattering matrix $S$.
- **Scattering Matrix (S-matrix)**: Uses a wavefunction matching technique (wavefront propagation) which is numerically stable for large systems.
- **Recursive Green's Functions (RGF)**: For computing local quantities like charge and spin densities.
- **Disorder Averaging**: Efficiently handles random disorder via supercell averaging or configuration sampling.

## Capabilities
- **Transport Properties**:
  - Spin-dependent conductance ($G_{\uparrow}, G_{\downarrow}$).
  - Transmission coeffecients $T(E, \mathbf{k}_{||})$.
  - Shot noise and Fano factor.
- **Spintronics**:
  - Spin-Transfer Torque (STT).
  - Spin-Orbit Torque (SOT).
  - Spin diffusion capability.
  - Magnetocrystalline anisotropy (if SOC included).
- **System Types**:
  - Magnetic Tunnel Junctions (MTJs).
  - Spin Valves (GMR/TMR).
  - Domain Walls and Skrymions.
  - Point Contacts.

## Key Strengths
- **Supercell Scalability**: Unlike CPA codes, TQT explicitly treats disorder in large lateral supercells, capturing effects like Anderson localization and diffusive transport.
- **Versatility**: Unified treatment of TB and MTO Hamiltonians.
- **Stability**: The S-matrix implementation is extremely stable for long systems where standard transfer matrix methods fail.
- **Non-Collinear Magnetism**: Full support for non-collinear spin textures (domain walls, spin spirals).

## Inputs & Outputs
- **Inputs**:
  - Hamiltonian files (TB or LMTO/EMTO format).
  - `input` config file: Geometry, energy range, k-points, disorder settings.
- **Outputs**:
  - `conductance.dat`: Energy/k-resolved transmission.
  - `density.dat`: Local density of states/charge.
  - `currents.dat`: Spin/Charge current distributions.

## Interfaces & Ecosystem
- **Upstream**:
  - **Crary** (LMTO code often used at Twente).
  - **EMTO**: Can map EMTO parameters to tight-binding forms.
- **Analysis**: Python processing tools provided in the repository.

## Performance Characteristics
- **Computational Cost**: $O(N)$ scaling with system length (thanks to RGF/S-matrix).
- **Parallelism**: MPI parallelization over energy points and transverse k-points ($k_x, k_y$). Efficient for high-throughput screening.
- **Memory**: Moderate; stores slice-by-slice matrices, avoiding full system Hamiltonian storage.

## Limitations & Known Constraints
- **Electrostatics**: Typically uses a "frozen potential" or simple self-consistency; possibly less rigorous Poisson solving than NEGF-DFT codes like TranSIESTA.
- **Basis**: Relies on localized orbital descriptions; no plane-wave support.

## Comparison with Other Codes
- **vs. Kwant**: Kwant is Python-based and very flexible for models; TQT is optimized Fortran 95/2003 with deeper support for specific *ab initio* MTO bases and spintronic observables.
- **vs. Smeagol**: Smeagol is DFT-NEGF; TQT is often used as a "post-DFT" transport solver on fitted or MTO Hamiltonians, allowing larger system sizes.

## Application Areas
- **MRAM**: Modeling tunnel magnetoresistance in Fe/MgO/Fe junctions.
- **Spin Logic**: Spin-orbit torque switching simulations.
- **Material Science**: Scattering by grain boundaries and interface roughness in metals (Cu interconnects).

## Community and Support
- **Development**: Developed at the University of Twente (Kelly, Starikov groups).
- **Source**: GitHub.

## Verification & Sources
- **Repository**: [https://github.com/TwenteQT/TwenteQuantumTransport](https://github.com/TwenteQT/TwenteQuantumTransport)
- **Primary Publication**: Found in repo README (Starikov et al.).
- **Verification status**: ✅ VERIFIED
  - Active repo with recent commits.
  - Well-documented functionality.
