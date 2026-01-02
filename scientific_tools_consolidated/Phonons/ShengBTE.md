# ShengBTE

## Official Resources
- Homepage: http://www.shengbte.org/
- Documentation: http://www.shengbte.org/index.php/usage
- Source Repository: http://www.shengbte.org/downloads
- License: GNU General Public License v3.0

## Overview
ShengBTE is a software package for computing lattice thermal conductivity of crystalline materials from first principles by solving the Boltzmann transport equation (BTE) for phonons. It uses harmonic and anharmonic interatomic force constants from DFT calculations.

**Scientific domain**: Thermal transport, thermoelectrics, phonon physics  
**Target user community**: Researchers studying thermal conductivity and phonon-mediated properties

## Theoretical Methods
- Boltzmann transport equation (BTE) for phonons
- Iterative and relaxation time approximation solutions
- Third-order anharmonic force constants
- Three-phonon scattering processes
- Normal and Umklapp processes
- Isotope scattering
- Boundary scattering
- Finite-size effects

## Capabilities (CRITICAL)
- Lattice thermal conductivity from first principles
- Temperature-dependent thermal conductivity
- Phonon relaxation times
- Phonon mean free paths
- Phonon group velocities
- Mode-resolved thermal conductivity contributions
- Spectral thermal conductivity
- Cumulative thermal conductivity
- Grain boundary scattering effects
- Nanostructure size effects
- Isotope disorder scattering
- Anisotropic thermal conductivity tensors

**Sources**: Official ShengBTE documentation, cited in 6/7 source lists

## Inputs & Outputs
- **Input formats**:
  - CONTROL file (calculation parameters)
  - FORCE_CONSTANTS_2ND (harmonic)
  - FORCE_CONSTANTS_3RD (anharmonic)
  - From thirdorder.py or other phonon codes
  
- **Output data types**:
  - BTE.kappa (thermal conductivity tensor)
  - BTE.KappaTensor (full output)
  - BTE.ReciprocalLattice
  - Mode-resolved properties
  - Cumulative thermal conductivity data

## Interfaces & Ecosystem
- **Force constant generation**:
  - thirdorder.py - automated 3rd order FC calculation
  - Phonopy for 2nd order
  - Interfaces with VASP, Quantum ESPRESSO, etc.
  
- **Post-processing**:
  - Python scripts for analysis
  - Plotting utilities
  - cumulative_kappa.py for analysis

## Limitations & Known Constraints
- **Requires force constants**: Needs 2nd and 3rd order from DFT
- **Computational cost**: 3rd order FC calculations expensive
- **Supercell size**: Large supercells needed for accurate 3rd order FC
- **Memory**: Can be intensive for large systems
- **Classical approximation**: Uses classical phonon occupations (fixable)
- **Three-phonon only**: Higher-order processes not included
- **Learning curve**: Requires understanding of phonon BTE
- **Platform**: Linux/Unix

## Verification & Sources
**Primary sources**:
1. Official website: http://www.shengbte.org/
2. Documentation: http://www.shengbte.org/index.php/usage
3. W. Li et al., Comput. Phys. Commun. 185, 1747 (2014) - ShengBTE paper
4. thirdorder.py companion code

**Secondary sources**:
1. ShengBTE tutorials and examples
2. Published thermal conductivity calculations
3. Thermoelectric material studies
4. Confirmed in 6/7 source lists (claude, g, gr, k, m, q)

**Confidence**: CONFIRMED - Appears in 6 of 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: ACCESSIBLE
- Source code: Available (download from website)
- Community support: Active (email, publications)
- Academic citations: >500
- Standard tool: Widely used for thermal conductivity
