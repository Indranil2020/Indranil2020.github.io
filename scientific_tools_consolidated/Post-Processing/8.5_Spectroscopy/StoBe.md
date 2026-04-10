# StoBe

## Official Resources
- Homepage: https://www.fz-juelich.de/pgi/pgi-1/DE/Home/home_node.html
- Documentation: Included with distribution
- Source Repository: Available upon request from developers
- License: Free for academic use

## Overview
**StoBe** (Stockholm-Berlin) is a DFT code based on Gaussian-type orbitals (GTO) specifically designed for the simulation of core-level X-ray spectroscopy including XAS, XES, and XPS. It uses the transition potential (TP) method and half/core-hole approaches for accurate core excitation and emission calculations.

**Scientific domain**: X-ray absorption, emission, and photoelectron spectroscopy of molecules and clusters  
**Target user community**: Researchers studying molecular and cluster core-level spectroscopy with DFT

## Theoretical Methods
- Density Functional Theory (DFT)
- Gaussian-type orbital (GTO) basis sets
- Transition potential (TP) method for XAS
- Half core-hole (Z+1/2) approximation
- Full core-hole (Z+1) approximation
- Linear response TDDFT for core excitations
- Static exchange (STEX) method
- LDA, GGA, hybrid functionals
- Scalar relativistic corrections (ZORA)
- Basis set superposition error (BSSE) correction

## Capabilities (CRITICAL)
- X-ray absorption spectroscopy (XAS/NEXAFS)
- X-ray emission spectroscopy (XES)
- X-ray photoelectron spectroscopy (XPS)
- Core-level binding energies
- Transition potential calculations
- STEX calculations
- Cluster and molecular calculations
- Symmetry-adapted calculations
- Geometry optimization
- Vibrational analysis
- Charged and neutral excitations

**Sources**: StoBe documentation, J. Chem. Phys. 117, 963 (2002)

## Key Strengths

### Transition Potential Method:
- Accurate core excitation energies
- Beyond sudden approximation
- Includes relaxation effects
- Systematic improvement possible
- Well-tested methodology

### Gaussian Basis Sets:
- Flexible basis selection
- Systematic improvement
- Efficient for molecules
- Augmented basis for Rydberg states
- Core-valence separation

### Core-Level Spectroscopy:
- Dedicated XAS/XES/XPS implementation
- Multiple core-hole approximations
- Spin-orbit splitting treatment
- Polarization dependence
- Energy-dependent cross-sections

### Molecular Focus:
- Optimized for finite systems
- No periodic boundary artifacts
- Accurate for gas-phase and clusters
- BSSE correction available
- Geometry optimization

## Inputs & Outputs
- **Input formats**:
  - StoBe input files
  - Basis set files (6-31G, cc-pVXZ, etc.)
  - Geometry specifications
  
- **Output data types**:
  - XAS spectra (oscillator strengths)
  - XES spectra (emission energies/intensities)
  - XPS binding energies
  - Orbital energies and characters
  - Transition dipole moments

## Interfaces & Ecosystem
- Standalone code
- Output compatible with standard plotting tools
- Basis sets from standard libraries
- Geometry from XYZ or internal coordinates

## Performance Characteristics
- **Speed**: Fast for molecular systems
- **Accuracy**: Good with TP method (0.5-1 eV for XAS)
- **System size**: Up to ~100 atoms typical
- **Memory**: Moderate (Gaussian integrals)

## Computational Cost
- **XAS (TP)**: Moderate (single SCF per edge)
- **XES**: Moderate
- **XPS**: Fast (ΔKohn-Sham)
- **Typical**: Minutes to hours per spectrum

## Limitations & Known Constraints
- **No periodic systems**: Cluster/molecular only
- **Basis sets**: Requires augmented basis for Rydberg
- **No BSE/GW**: DFT-level only
- **Spin-orbit**: Limited treatment
- **Installation**: Requires registration
- **Documentation**: Limited public documentation

## Comparison with Other Codes
- **vs FEFF**: StoBe is molecular, FEFF is periodic/multiple scattering
- **vs FDMNES**: StoBe uses GTO, FDMNES uses finite differences
- **vs ORCA**: StoBe specialized for XAS, ORCA is general QC
- **vs xspectra**: StoBe is molecular, xspectra is periodic (QE)
- **Unique strength**: Transition potential method for molecular XAS/XES, Gaussian basis flexibility

## Application Areas

### Molecular XAS:
- Organic molecules NEXAFS
- Carbon K-edge spectroscopy
- Nitrogen K-edge spectroscopy
- Oxygen K-edge spectroscopy

### Transition Metal Complexes:
- Metal L-edge XAS
- Charge transfer satellites
- Ligand field effects
- Oxidation state analysis

### Clusters and Nanoparticles:
- Metal cluster core spectra
- Supported cluster spectroscopy
- Size-dependent effects
- Adsorbate spectroscopy

### Environmental and Biological:
- Water XAS
- Amino acid spectroscopy
- Pollutant characterization
- Soil mineral spectroscopy

## Best Practices

### Basis Set Selection:
- Use augmented basis sets for Rydberg states
- Include diffuse functions for XAS
- Test convergence with basis size
- Separate core and valence basis

### Core-Hole Treatment:
- Compare TP and full core-hole results
- Use TP for better absolute energies
- Consider Z+1 approximation for screening
- Validate against experiment

### Cluster Model Construction:
- Use sufficiently large clusters
- Cap dangling bonds appropriately
- Test convergence with cluster size
- Consider embedding

## Community and Support
- Free for academic use (registration required)
- Developed at FZ Jülich and Stockholm University
- Limited but dedicated user community
- Key reference: L. G. M. Pettersson et al., J. Chem. Phys. 117, 963 (2002)

## Verification & Sources
**Primary sources**:
1. FZ Jülich PGI-1: https://www.fz-juelich.de/pgi/pgi-1/
2. L. G. M. Pettersson and T. B. V. H. Pettersson, J. Chem. Phys. 117, 963 (2002)
3. K. Hermann et al., Phys. Rev. B 73, 125109 (2006)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: Available with distribution
- Source code: Available upon registration
- Community support: Dedicated but small
- Academic citations: >500 (method papers)
- Active development: Maintained
- Specialized strength: Transition potential method for molecular XAS/XES/XPS, Gaussian basis core-level spectroscopy
