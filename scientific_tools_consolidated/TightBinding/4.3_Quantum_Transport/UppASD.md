# UppASD

## Official Resources
- **Homepage**: https://uppasd.github.io/
- **Repository**: https://github.com/UppASD/UppASD
- **License**: GPL-3.0

## Overview
**UppASD** (Uppsala Atomistic Spin Dynamics) is a premier software package for performing **atomistic spin dynamics (ASD)** simulations. It models the time evolution of magnetic moments in materials at **finite temperature** by solving the stochastic Landau-Lifshitz-Gilbert (LLG) equation. It effectively treats the magnetic system as a lattice of classical spins interacting via Heisenberg exchange, DMI, and anisotropy, typically parameterized from first-principles (DFT) calculations.

**Scientific domain**: Magnetism, Spintronics, Thermodynamics
**Target user community**: Researchers bridging DFT and finite-temperature magnetism

## Theoretical Methods
- **Stochastic LLG**: Integrates the equation of motion $\frac{d\mathbf{m}}{dt} = - \gamma \mathbf{m} \times (\mathbf{B}_{eff} + \mathbf{b}_{fl}) + \alpha \mathbf{m} \times \frac{d\mathbf{m}}{dt}$.
- **Monte Carlo**: Metropolis MC for equilibrium thermodynamic properties.
- **Hamiltonian**:
  - Heisenberg Exchange ($J_{ij}$).
  - Dzyaloshinskii-Moriya Interaction (DMI).
  - Magnetic Anisotropy ($K$).
  - Dipole-Dipole interactions.
  - Higher-order (biquadratic, 4-spin) terms.
- **Spin-Lattice Dynamics**: Coupled simulation of spin and lattice (phonon) degrees of freedom (in recent versions).

## Capabilities
- **Thermodynamics**:
  - Curie/Néel temperatures ($T_C, T_N$).
  - Magnetization vs Temperature curves $M(T)$.
  - Magnetic phase diagrams.
- **Dynamics**:
  - Spin wave dispersion and lifetimes (via Dynamic Structure Factor $S(\mathbf{q}, \omega)$).
  - Ultrafast demagnetization dynamics.
  - Domain wall motion and skyrmion diffusion.
- **Correlation**:
  - Space-time correlation functions.

## Key Strengths
- **Material Specificity**: Designed to work with parameters mapped directly from DFT codes (like RSPt, KKR, VASP), making it highly predictive for real materials (Fe, Co, NiO, etc.).
- **Magnonics**: Excellent for calculating magnon spectra and lifetimes to compare with neutron scattering (INS) experiments.
- **Efficiency**: Highly optimized Fortran code capable of simulating millions of spins.

## Inputs & Outputs
- **Inputs**:
  - `j` content files: Lists of exchange parameters ($J_{ij}$).
  - `inp`: Main control file (lattice, temperature, steps).
- **Outputs**:
  - Moment files (`moment.out`).
  - Correlation functions (`sqw.out`).
  - Thermodynamics (`averages.out`).

## Interfaces & Ecosystem
- **Upstream**: Strong links with DFT codes like **RSPt** (Uppsala) and **SPR-KKR** (Munich) for extracting $J_{ij}$.
- **Visualization**: Output compatible with ParaView/OVITO.

## Performance Characteristics
- **Computational Cost**: $O(N)$ for short-range interactions.
- **Parallelism**: Hybrid MPI/OpenMP parallelization; GPU (CUDA) support available for massive systems.

## Comparison with Other Codes
- **vs. Spirit**: UppASD traditionally focuses more on thermodynamic averages and spectra ($S(q, \omega)$); Spirit excels at finding energy barriers (GNEB) and interactive dynamics.
- **vs. Vampire**: Very similar scope; UppASD has historically strong connections to the ab initio community in Europe (Psi-k).

## Community and Support
- **Development**: Uppsala University (Olle Eriksson, B. Skubic) and collaborators.
- **Source**: GitHub.

## Verification & Sources
- **Website**: [https://uppasd.github.io/](https://uppasd.github.io/)
- **Primary Publication**: B. Skubic et al., J. Phys.: Condens. Matter 20, 315203 (2008).
- **Verification status**: ✅ VERIFIED
  - A staple code in the computational magnetism community.
