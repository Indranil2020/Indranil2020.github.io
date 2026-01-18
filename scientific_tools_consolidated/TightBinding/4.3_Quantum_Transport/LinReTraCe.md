# LinReTraCe

## Official Resources
- **Repository**: https://github.com/linretrace/linretrace
- **Documentation**: https://github.com/linretrace/linretrace/wiki
- **License**: GPL-3.0

## Overview
**LinReTraCe** (Linear Response Transport Centre) is a package for the simulation of linear response transport properties in solids. It employs a **semi-analytical evaluation of Kubo formulas** (such as the Kubo-Bastin formula) to determine resistivities and coefficients for Hall, Seebeck, and Nernst effects. Uniquely, it is designed to capture the effects of electronic correlations and finite quasiparticle lifetimes by accepting **self-energy** inputs from methods like DMFT (Dynamical Mean-Field Theory) or many-body perturbation theory.

**Scientific domain**: Correlated Electrons, Topological Transport, Thermoelectrics
**Target user community**: Theorists working on strongly correlated systems and topological materials

## Theoretical Methods
- **Kubo Formula**: Evaluation of the current-current correlation function to obtain conductivity tensors.
- **Semi-Analytical Integration**: Performs analytical integration over Green's functions products where possible, allowing for extremely high precision.
- **Self-Energy ($\Sigma(\omega)$)**: Input-driven scattering rates, allowing modeling of impurity scattering, electron-electron correlation (DMFT), or electron-phonon coupling.
- **Tetrahedron Method**: Adaptive integration over the Brillouin zone for Fermi surface properties.

## Capabilities
- **Coefficients**:
  - DC Conductivity/Resistivity ($\sigma_{xx}, \rho_{xx}$).
  - Hall Coefficient ($R_H$) and Anomalous Hall Conductivity ($\sigma_{xy}$).
  - Seebeck Coefficient ($S$).
  - Nernst Signal ($\nu$).
  - Optical Conductivity $\sigma(\omega)$.
- **Analysis**:
  - Band-resolved contributions.
  - Fermi surface visualization (spectral weight).

## Key Strengths
- **Convergence**: The semi-analytical approach allows for the use of extremely dense k-grids ($1000^3$), which is often necessary to resolve fine structures like the Nernst signal in topological materials.
- **Correlations**: Seemlessly integrates with DMFT outputs, making it a standard tool for transport in heavy fermions and cuprates.
- **Versatility**: Interfaces with Wannier90, making it compatible with virtually all DFT codes.

## Inputs & Outputs
- **Inputs**:
  - Hamiltonian: `wannier90_hr.dat` or internal models.
  - Self-Energy: Text files containing $\Sigma(\omega)$.
  - Chemical potential and Temperature range.
- **Outputs**:
  - Text files with temperature/doping dependent coefficients.

## Interfaces & Ecosystem
- **Wannier90**: Primary interface for realistic materials.
- **DFT Codes**: VASP, Wien2k, FPLO (via Wannier90).
- **DMFT Codes**: Can read self-energies compatible with standard formats.

## Performance Characteristics
- **Speed**: Highly optimized C++ core with MPI/OpenMP.
- **Accuracy**: Superior to standard "constant relaxation time" Boltzmann codes for anomalous/topological quantities.

## Comparison with Other Codes
- **vs. BoltzTrap**: BoltzTrap uses Boltzmann theory (semi-classical); LinReTraCe uses Kubo (quantum mechanical) and handles self-energies/scattering explicitly (beyond constant $\tau$).
- **vs. Wannier90 (Berry)**: Wannier90 considers intrinsic Berry curvature terms; LinReTraCe includes extrinsic scattering contributions (vertex corrections) more naturally via the Green's function.

## Application Areas
- **Strange Metals**: Transport in non-Fermi liquids.
- **Topological Semimetals**: Weyl/Dirac points signatures in Hall/Nernst effects.
- **Thermoelectrics**: Optimization of power factors in complex oxides.

## Community and Support
- **Development**: TU Dresden and University of Augsburg.
- **Source**: GitHub.

## Verification & Sources
- **Repository**: [https://github.com/linretrace/linretrace](https://github.com/linretrace/linretrace)
- **Primary Publication**: Refer to repository README for latest citation (e.g., *Computer Physics Communications*).
- **Verification status**: ✅ VERIFIED
  - Active academic project.
