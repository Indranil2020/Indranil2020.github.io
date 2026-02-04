# xspectra

## Official Resources
- Homepage: https://www.quantum-espresso.org/
- Documentation: https://www.quantum-espresso.org/Doc/xspectra_user_guide/
- Source Repository: https://gitlab.com/QEF/q-e (Part of Quantum ESPRESSO)
- License: GNU General Public License v2.0

## Overview
XSpectra is a code for calculating X-ray Absorption Spectra (XAS) at the K-edge (and L-edge) using the Projector Augmented Wave (PAW) method or Pseudopotentials. It is part of the Quantum ESPRESSO distribution (`xspectra.x`). It avoids the explicit calculation of empty states by using the Lanczos recursion algorithm to compute the continued fraction representation of the Green's function.

**Scientific domain**: X-ray spectroscopy, XANES, core-level excitations  
**Target user community**: Spectroscopists, Quantum ESPRESSO users

## Theoretical Methods
- X-ray Absorption Near Edge Structure (XANES)
- Dipole and Quadrupole approximation
- Lanczos recursion method (continued fraction)
- Projector Augmented Wave (PAW) reconstruction of all-electron wavefunction
- Core-hole effects (via supercells with core-hole pseudopotentials)

## Capabilities (CRITICAL)
- Calculation of K-edge and L2,3-edge XAS
- Linear and circular dichroism
- Efficient calculation for large systems (no empty states needed)
- Dependence on polarization
- Core-hole treatment (FCH/XCH approximations)
- Interface with `pw.x` charge density

**Sources**: XSpectra documentation, Phys. Rev. B 80, 035102 (2009)

## Inputs & Outputs
- **Input formats**: `xspectra.in` (namelist), `prefix.save` (from pw.x)
- **Output data types**: `xspectra.dat` (energy vs cross-section)

## Interfaces & Ecosystem
- **Quantum ESPRESSO**: Fully integrated
- **GIPAW**: Uses GIPAW reconstruction
- **Visualization**: Output is plain text for plotting

## Workflow and Usage
1. Perform SCF calculation with `pw.x` (ground state or with core-hole).
2. Create `xspectra.in`: Define edge, absorbing atom, Lanczos parameters.
3. Run `xspectra.x`.
4. Plot the resulting spectrum.

## Performance Characteristics
- Highly efficient due to Lanczos algorithm (scales linearly with N)
- Memory efficient compared to sum-over-states methods

## Application Areas
- Structure determination
- Oxidation state analysis
- Surface adsorption geometry
- High-pressure phases

## Community and Support
- Part of Quantum ESPRESSO community
- Active mailing list
- Standard tool for XANES

## Verification & Sources
**Primary sources**:
1. Homepage: https://www.quantum-espresso.org/
2. Publication: O. Bunau and M. Calandra, Phys. Rev. B 80, 035102 (2009)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: AVAILABLE
- Source: OPEN (GPL, part of QE)
- Development: ACTIVE (QE Foundation)
- Applications: XANES, PAW, Lanczos, core-hole
