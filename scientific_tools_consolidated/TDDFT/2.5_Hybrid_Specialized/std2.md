# std2 (Simplifed TD-DFT)

## Official Resources
- Homepage: https://github.com/grimme-lab/stda
- Source Repository: https://github.com/grimme-lab/stda
- License: GNU General Public License v3.0

## Overview
std2 is a software package developed by the Grimme group for performing simplified Time-Dependent Density Functional Theory (sTD-DFT) and simplified Tamm-Dancoff Approximation (sTDA) calculations. These methods provide a highly efficient approximation to full TDDFT, allowing for the computation of ultra-fast UV-Vis absorption and electronic circular dichroism (ECD) spectra for very large molecular systems (thousands of atoms).

**Scientific domain**: Large-scale excited states, high-throughput screening, UV-Vis/ECD spectroscopy
**Target user community**: Researchers screening large molecules, supramolecular chemists

## Theoretical Methods
- Simplified TDA (sTDA)
- Simplified TD-DFT (sTD-DFT)
- Simplified TDA-xTB (sTDA-xTB)
- Tight-binding approximations
- Multipole approximations for integrals
- Point-charge approximations
- Range-separated hybrid functionality

## Capabilities (CRITICAL)
- Excitation energies (singlet and triplet)
- Oscillator strengths (UV-Vis)
- Rotatory strengths (ECD)
- High-throughput spectral calculation
- Systems with 1000-5000 atoms
- Full spectral range (Valence/Rydberg)
- Analysis of excitons

**Sources**: GitHub repository, Grimme group publications

## Key Strengths

### Extreme Speed:
- Orders of magnitude faster than standard TDDFT
- Seconds/minutes for large molecules
- Integral approximations avoid N^4 scaling

### Large Scale:
- Routinely handles >1000 atoms
- Supramolecular complexes
- Protein fragments
- Nanoclusters

### Accuracy:
- Errors ~0.2-0.4 eV (comparable to range-separated hybrids)
- Calibrated against high-level methods
- Reliable spectral shapes

## Inputs & Outputs
- **Input formats**:
  - `molden` file (from ORCA, Turbomole, Gaussian, etc.)
  - `xtb` output files (for sTDA-xTB)
  - Control input file
  
- **Output data types**:
  - Excitation energies
  - Oscillator/Rotatory strengths
  - Simulated spectra (broadened)
  - State character analysis

## Interfaces & Ecosystem
- **QC Codes (for wavefunction)**: ORCA, TURBOMOLE, Gaussian, Q-Chem, PSI4 (via Molden)
- **xTB**: Seamless integration for sTDA-xTB
- **Language**: Fortran
- **Binaries**: Static binaries available

## Advanced Features

### xTB Integration:
- Can run purely on xTB wfn (sTDA-xTB)
- No DFT calculation required
- Extremely fast workflow

### Solvent Effects:
- Implicit solvation models accessible via interfaced code
- ALPB (in xTB)

## Performance Characteristics
- **Speed**: Ultra-fast (seconds)
- **Accuracy**: Qualitative to semi-quantitative
- **System size**: Up to 5000+ atoms
- **Memory**: Efficient storage of approximated integrals

## Computational Cost
- **Wavefunction**: Requires ground state (DFT or xTB)
- **Excited State**: Negligible compared to DFT
- **Scaling**: N^2 or better with approximations

## Limitations & Known Constraints
- **Approximation**: Not ab initio TDDFT
- **Wavefunction dependency**: Quality depends on input orbitals
- **Charge Transfer**: Corrected sTDA can handle it, but check
- **Rydberg states**: Can be limited by basis set
- **Input**: Requires Molden file (except xTB mode)

## Comparison with Other Codes
- **vs Full TDDFT**: std2 is approx. 100-1000x faster, less rigorous
- **vs ZINDO**: std2 generally more robust and accurate
- **vs DFTB**: Similar niche, std2 uses DFT orbitals
- **Unique strength**: Unmatched speed for realistic TDDFT-quality spectra of huge variants

## Application Areas
- **Screening**: Calculating spectra for thousands of conformers
- **Supramolecular**: Host-guest complexes
- **Bio-organic**: Large chromophores in proteins
- **Materials**: Optical properties of large aggregates

## Best Practices
- **Input Orbitals**: Use robust DFT functionals (e.g. wB97X-D)
- **Basis Set**: def2-SVP/TZVP usually sufficient
- **TDA vs TDDFT**: sTDA usually robust, sTD-DFT includes de-excitation
- **Verification**: Check a small model with full TDDFT

## Community and Support
- Open-source GPL v3
- Grimme group support
- Binaries provided
- Used in `xtb` ecosystem

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/grimme-lab/stda
2. S. Grimme, J. Chem. Phys. 138, 244104 (2013)
3. C. Bannwarth, S. Grimme, Comput. Theor. Chem. 1040, 45 (2014)

**Confidence**: VERIFIED - Grimme group official code

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Source code: OPEN (GPL v3)
- Method: sTDA / sTD-DFT (Widely cited)
- Specialized strength: Ultra-fast excited states for massive systems
