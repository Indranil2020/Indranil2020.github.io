# Libwfa

## Official Resources
- Homepage: https://github.com/libwfa/libwfa
- Documentation: https://libwfa.github.io/libwfa/
- Source Repository: https://github.com/libwfa/libwfa
- License: BSD 3-Clause License

## Overview
Libwfa is an open-source C++ library for wavefunction analysis of electronic excitations. It implements various methods for analyzing excited state calculations (e.g., from TD-DFT, ADC, CC, EOM-CC) by computing state difference density matrices, natural transition orbitals (NTOs), and charge-transfer descriptors. It provides tools to visualize and quantify the character of electronic transitions.

**Scientific domain**: Wavefunction analysis, excited states, charge transfer  
**Target user community**: Quantum chemists, photochemists, materials scientists

## Theoretical Methods
- Natural Transition Orbitals (NTOs)
- State Difference Density Matrices (SDDM)
- Detachment/Attachment Density Matrices
- Charge Transfer Numbers (q_CT)
- Exciton size and electron-hole correlation
- Mulliken and Lowdin population analysis for excitations

## Capabilities (CRITICAL)
- Automated analysis of excited states
- Quantitative descriptors for charge transfer (CT) character
- Visualization of electron-hole pairs via NTOs
- Interface with quantum chemistry codes (Q-Chem, MOLCAS, Orca, etc.)
- Calculation of exciton binding energies (qualitative)
- Decomposition of excitation energy

**Sources**: Libwfa documentation, J. Comput. Chem. 37, 1632 (2016)

## Inputs & Outputs
- **Input formats**: Code-specific wavefunction/density data (HDF5 or native formats)
- **Output data types**: Analysis log, cube files for NTOs/densities, CT descriptors

## Interfaces & Ecosystem
- **Q-Chem**: Integrated directly
- **OpenMolcas**: Integrated directly
- **Orca**: Can be used via interfaces
- **Molden/Cube**: Visualization output

## Workflow and Usage
1. Perform excited state calculation (e.g., TD-DFT).
2. Generate required density matrices (Transition, Difference).
3. Run Libwfa analysis (often built-in to the host code).
4. Visualize NTOs using VMD or Molden.

## Performance Characteristics
- Highly efficient linear algebra operations
- Negligible cost compared to the excited state calculation itself

## Application Areas
- Photovoltaic materials (charge separation)
- OLEDs (TADF materials)
- Photocatalysis
- Biological light harvesting

## Community and Support
- Open-source (BSD)
- Developed by Plasser, Wormit, and Dreuw groups
- Active development and integration

## Verification & Sources
**Primary sources**:
1. Homepage: https://github.com/libwfa/libwfa
2. Publication: F. Plasser, M. Wormit, A. Dreuw, J. Chem. Phys. 141, 024106 (2014)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE (GitHub)
- Documentation: AVAILABLE
- Source: OPEN (GitHub)
- Development: ACTIVE
- Applications: NTOs, excited state analysis, charge transfer
