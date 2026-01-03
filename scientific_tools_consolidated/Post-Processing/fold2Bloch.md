# fold2Bloch

## Official Resources
- Homepage: https://github.com/rubel75/fold2Bloch-VASP
- Documentation: https://github.com/rubel75/fold2Bloch-VASP/wiki
- Source Repository: https://github.com/rubel75/fold2Bloch-VASP
- License: MIT License

## Overview
fold2Bloch is a code for unfolding the band structure of supercells obtained from first-principles calculations. Specifically designed for VASP (and adaptable to other plane-wave codes), it recovers the Bloch character of electronic eigenstates in the primitive Brillouin zone. This allows for the direct comparison of supercell calculations (used for defects, alloys, or magnetic structures) with the band structures of primitive cells and experimental ARPES data.

**Scientific domain**: Band unfolding, supercell calculations, electronic structure  
**Target user community**: VASP users, defect researchers, alloy physicists

## Theoretical Methods
- Band unfolding algorithm
- Spectral weight calculation
- Projection of supercell wavefunctions onto primitive cell k-points
- Recovery of Bloch character ($P_{Km}$)
- Analysis of symmetry breaking effects

## Capabilities (CRITICAL)
- Unfolding of electronic bands from supercells
- Calculation of spectral weights for primitive k-points
- Identification of "ghost bands" vs real bands
- Support for non-spin-polarized and spin-polarized calculations
- Handling of large supercells
- Output formatting for plotting (gnuplot, etc.)

**Sources**: fold2Bloch documentation, Phys. Rev. B 81, 121201(R) (2010)

## Inputs & Outputs
- **Input formats**: WAVECAR (VASP wavefunction), fold2Bloch input file (defining unfolding transformation)
- **Output data types**: Spectral weights (typically `*.f2b` files), formatted data for plotting band structures

## Interfaces & Ecosystem
- **VASP**: Primary interface (requires WAVECAR)
- **Wien2k**: A version exists for Wien2k (`fold2Bloch-Wien2k`)
- **Plotting**: Output compatible with standard plotters

## Workflow and Usage
1. Perform VASP calculation on supercell (generate WAVECAR).
2. Determine transformation matrix between supercell and primitive cell.
3. Run fold2Bloch: `fold2Bloch WAVECAR output_name`
4. Process output to generate band structure plot (energy vs k-path weighted by character).

## Performance Characteristics
- Fast post-processing step
- Memory requirement depends on WAVECAR size
- Efficient for large supercells compared to some other methods

## Application Areas
- Defect levels in semiconductors
- Band gap renormalization in alloys
- Interface states
- Disordered systems modeled by supercells

## Community and Support
- Open-source (MIT)
- GitHub repository
- Developed by O. Rubel group (McMaster University)

## Verification & Sources
**Primary sources**:
1. GitHub (VASP): https://github.com/rubel75/fold2Bloch-VASP
2. GitHub (Wien2k): https://github.com/rubel75/fold2Bloch
3. Publication: V. Popescu and A. Zunger, Phys. Rev. B 85, 085201 (2012) (Methodology ref)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: AVAILABLE
- Source: OPEN (GitHub)
- Development: ACTIVE (Rubel Group)
- Applications: Band unfolding, VASP/Wien2k, supercells
