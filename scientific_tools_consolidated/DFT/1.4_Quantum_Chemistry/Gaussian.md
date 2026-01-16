# Gaussian

## Official Resources
- Homepage: https://gaussian.com/
- Documentation: https://gaussian.com/man/
- Source Repository: Proprietary (commercial license)
- License: Commercial license required

## Overview
Gaussian is the most widely-used electronic structure program in chemistry. Originally developed by John Pople (Nobel Prize 1998) and now maintained by Gaussian, Inc., it provides a comprehensive suite of methods from semi-empirical to high-level coupled cluster with extensive automation and user-friendly interface. The current version is Gaussian 16, with GaussView 6 as the companion visualization tool.

**Scientific domain**: Computational chemistry, drug design, materials chemistry, spectroscopy  
**Target user community**: Chemists across academia and industry; pharmaceutical and materials companies

## Theoretical Methods
- Hartree-Fock (RHF, UHF, ROHF)
- Semi-empirical methods (AM1, PM3, PM6, PM7)
- Density Functional Theory (DFT)
  - LDA, GGA, meta-GGA, hybrid, double-hybrid functionals
  - B3LYP, PBE, TPSS, APFD, MN15, MN15L, ωB97X-D
- Møller-Plesset (MP2, MP3, MP4, MP5)
- Coupled Cluster (CCSD, CCSD(T))
- Configuration Interaction (CI, QCISD, QCISD(T))
- Complete Active Space (CASSCF, CASMP2)
- Time-Dependent DFT (TDDFT)
- EOM-CCSD for excited states
- Composite methods (CBS-QB3, CBS-APNO, G4, W1, W2)
- Solvation models (PCM, SMD, CPCM)
- Dispersion corrections (GD2, GD3, GD3BJ)
- ONIOM for QM/MM and multi-layer calculations
- Relativistic ECPs (Stuttgart-Dresden, Ahlrichs)

## Capabilities (CRITICAL)
- Ground-state electronic structure
- Geometry optimization and conformational searches
- Transition state searches (QST2, QST3, Berny)
- IRC (Intrinsic Reaction Coordinate) calculations
- Vibrational frequencies and thermochemistry
- Anharmonic vibrational analysis
- Excited states (TDDFT, CIS, EOM-CCSD)
- Spectroscopic properties:
  - UV-Vis and fluorescence
  - IR and Raman spectra (including anharmonic)
  - NMR chemical shifts and coupling constants
  - EPR/ESR parameters
  - VCD (vibrational circular dichroism)
  - ROA (Raman optical activity)
  - Electronic circular dichroism (ECD)
  - CPL (circularly polarized luminescence)
  - Optical rotatory dispersion (ORD)
  - Resonance Raman
- Molecular properties (dipole, polarizability, hyperpolarizability)
- Potential energy surface scans
- ONIOM multi-layer QM/MM calculations
- Direct dynamics (ADMP, BOMD)
- Excitation Energy Transfer (EET)
- Solvation effects with analytical derivatives
- GPU acceleration (NVIDIA K40/K80 for HF/DFT)

**Sources**: Official Gaussian documentation, cited in 7/7 source lists

## Key Strengths

### Industry Standard:
- Most widely cited quantum chemistry code (>100,000 citations)
- De facto standard in pharmaceutical industry
- Established foundation for computational chemistry
- Extensive validation and benchmarking

### Automation:
- Extensive automated methods
- User-friendly input syntax
- Robust optimization algorithms
- Composite methods for thermochemistry
- Automated conformational analysis

### Comprehensive Methods:
- All major quantum chemical methods
- Extensive spectroscopy predictions
- Wide range of DFT functionals
- High-level correlation methods

### GaussView Integration:
- Intuitive graphical interface
- Visualization of results
- Input file builder
- Spectrum simulation

## Inputs & Outputs
- **Input formats**:
  - Route section with keywords
  - Z-matrix or Cartesian coordinates
  - Checkpoint files for restart
  - Simple, human-readable format
  
- **Output data types**:
  - Formatted output files (.log)
  - Checkpoint files (.chk)
  - Formatted checkpoint files (.fchk)
  - Cube files for densities and orbitals
  - Archive entries

## Interfaces & Ecosystem
- **Visualization**:
  - GaussView 6 - integrated GUI
  - Molden, Avogadro, Chemcraft compatible
  
- **Workflow integration**:
  - Widely supported by workflow tools
  - Python wrappers (cclib, GaussianWrangler)
  - ASE interface
  
- **Analysis tools**:
  - formchk - checkpoint file formatting
  - cubegen - cube file generation
  - freqchk - frequency analysis
  - c8616 - linking utilities


## Workflow and Usage

### Input Format:
Gaussian uses a route section to define the method and basis set, followed by the molecule specification.

```
#P B3LYP/6-31G(d) Opt Freq

Water Optimization and Frequency

0 1
 O
 H 1 0.96
 H 1 0.96 2 104.5

```

### Running Gaussian:
```bash
# Standard execution
g16 input.com
# Output goes to input.log by default
```

### Common Tasks:
- **Opt**: Geometry optimization
- **Freq**: Frequency analysis
- **SP**: Single point energy (default)
- **TD**: Excited states
- **SCRF**: Solvation models

## Advanced Features

### ONIOM (QM/MM):
- "Our own N-layered Integrated molecular Orbital and Molecular mechanics"
- Multi-layer calculations (High/Medium/Low)
- Treats large systems effectively
- Electronic embedding
- Automatic topology handling

### Composite Methods:
- High-accuracy thermochemistry
- CBS-QB3, G3, G4, W1
- Automated multi-step protocols
- Extrapolation to basis set limit
- Chemical accuracy (1 kcal/mol)

### Solvent Effects:
- SCRF (Self-Consistent Reaction Field)
- Polarizable Continuum Model (PCM)
- SMD (Solvation Model based on Density)
- State-specific solvation for excited states
- Non-equilibrium solvation

### Spectroscopic Prediction:
- VCD (Vibrational Circular Dichroism)
- ROA (Raman Optical Activity)
- NMR spin-spin coupling
- Anharmonic vibrational analysis
- Franck-Condon analysis

### Automated Transition State Search:
- QST2/QST3 (Synchronous Transit-Guided Quasi-Newton)
- Requires reactants and products (and TS guess for QST3)
- GDIIS algorithm
- Eigenvector-following

## Performance Characteristics
- **Speed**: Efficient integrals and SCF convergence
- **Scalability**: Shared-memory parallelization (OpenMP) efficient; Linda for distributed is limited compared to MPI codes like NWChem
- **GPU support**: Available for specific modules (DFT gradients/frequencies)
- **Memory**: Usage defined in input (`%Mem=NGB`), critical for performance
- **Disk I/O**: Heavy use of Read-Write files (.rwf)

## Computational Cost
- **DFT**: Efficient for medium systems (up to ~500 atoms)
- **MP2**: Moderate cost, O(N^5)
- **CCSD(T)**: Very expensive, O(N^7)
- **Composite Methods**: Expensive but automated
- **Frequencies**: Expensive (requires second derivatives)

## Comparison with Other Codes
- **vs ORCA**: Gaussian is commercial, better GUI (GaussView); ORCA is free for academia, better coupled cluster performance.
- **vs GAMESS**: Gaussian has more automation/composite methods; GAMESS is free/open-source.
- **vs NWChem**: Gaussian easier to use for small molecules; NWChem superior for massive parallelism.
- **vs Q-Chem**: Similar market; Q-Chem arguably stronger in recent DFT functionals and excited states.
- **Unique strength**: Ease of use, automation, huge user base, industry standard status, ONIOM method.

## Best Practices

### Input Configuration:
- Always check spin multiplicity
- Use `%NProcShared` to set processors
- Set `%Mem` appropriately (avoid swapping)
- Use `Opt=CalcFC` for transition states

### Basis Sets:
- Use Pople sets (6-31G*) for routine work
- Use Correlation Consistent (cc-pVTZ) for high accuracy
- Use Diffuse functions (+) for anions/excited states

### Convergence Issues:
- Use `SCF=XQC` for difficult cases
- Use `Opt=GDIIS` for shallow potentials
- Check wavefunction stability (`Stable=Opt`)

## Community and Support
- **Support**: Commercial support from Gaussian, Inc.
- **Resources**: "Exploring Chemistry with Electronic Structure Methods" (The "Gaussian Bible")
- **White Papers**: Technical details on website
- **Workshops**: Official training sessions
- **User Base**: Largest in the field, abundant online examples

## Application Areas

### Drug Discovery:
- Ligand-protein binding energies
- Conformational analysis
- Reactivity predictions
- ADMET property prediction

### Reaction Mechanisms:
- Transition state characterization
- Reaction pathways (IRC)
- Activation energies
- Thermochemistry with composite methods

### Spectroscopy:
- UV-Vis spectra
- IR/Raman for identification
- NMR chemical shift prediction
- Chiroptical properties (VCD, ECD)

## Limitations & Known Constraints
- **Commercial license**: Expensive; significant license fees
- **Closed source**: No source code access or modification
- **Molecular focus**: Not optimized for extended/periodic systems
- **Periodic systems**: Limited support
- **System size**: Practical limits ~500-1000 atoms for DFT
- **Parallelization**: Efficient but proprietary implementation
- **Platform**: Linux, macOS, Windows (commercial binaries)
- **License management**: Can be restrictive (node-locked, floating)
- **Integral limits**: Maximum atoms and basis functions in integral program
- **Updates**: Periodic major releases (not continuous)

## Verification & Sources
**Primary sources**:
1. Official website: https://gaussian.com/
2. Manual: https://gaussian.com/man/
3. M. J. Frisch et al., Gaussian 16, Revision C.01, Gaussian, Inc., Wallingford CT, 2016
4. Gaussian White Papers and Technical Notes
5. J. Pople citation (Nobel Prize 1998)

**Secondary sources**:
1. Gaussian tutorials and documentation
2. Published applications across all chemistry
3. Textbook references (standard in computational chemistry courses)
4. Confirmed in 7/7 source lists (claude, g, gr, k, m, q, z)

**Confidence**: CONFIRMED - Appears in all 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE and ACCESSIBLE
- License: Commercial (verified)
- Community support: Extensive (support, GaussView, tutorials)
- Academic citations: >100,000 (most cited quantum chemistry code)
- Industry standard: Dominant in pharmaceutical and materials industry
- Specialized strength: Comprehensive methods, automation, spectroscopy, industry standard
