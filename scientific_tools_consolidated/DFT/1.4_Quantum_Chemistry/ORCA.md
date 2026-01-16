# ORCA

## Official Resources
- Homepage: https://orcaforum.kofo.mpg.de/
- Documentation: https://www.faccts.de/docs/orca/current/
- Download: https://orcaforum.kofo.mpg.de/app.php/portal
- License: Free for academic use (registration required); commercial via FACCTs GmbH

## Overview
ORCA is a modern, general-purpose quantum chemistry program package featuring extensive capabilities in molecular electronic structure calculations. Developed by Frank Neese and coworkers at the Max Planck Institute für Kohlenforschung, ORCA is known for its user-friendly input, comprehensive methods (from semi-empirical to high-level coupled cluster), excellent spectroscopic property calculations, and efficient parallelization. Version 6.0 (released July 2024) represents a near-complete rewrite with major performance improvements.

**Scientific domain**: Quantum chemistry, molecular electronic structure, spectroscopy  
**Target user community**: Computational chemists, spectroscopists, transition metal researchers

## Theoretical Methods
- Hartree-Fock (RHF, UHF, ROHF)
- Semi-empirical methods (AM1, PM3, PM6, PM7, GFN-xTB)
- Density Functional Theory (LDA, GGA, meta-GGA, hybrid, double-hybrid)
  - Modern functionals: r²SCAN, ωB97M-V, PBE0, B3LYP, TPSS
  - Double-hybrid: DSDPBEP86, PBE0DH, PBEQIDH
- Range-separated functionals (CAM-B3LYP, ωB97X)
- Dispersion corrections (DFT-D3, DFT-D4)
- RIJCOSX acceleration for HF and hybrid-DFT
- Møller-Plesset perturbation theory (MP2, SCS-MP2, RI-MP2)
- Coupled Cluster (CCSD, CCSD(T), DLPNO-CCSD(T))
- Multi-reference methods (CASSCF, NEVPT2, MRCI, MRCPA)
- Time-Dependent DFT (TDDFT) with analytical gradients
- EOM-CCSD for excited states with analytical gradients
- Solvation models (CPCM, SMD with analytical Hessian)
- Relativistic methods (ZORA, DKH, X2C)
- Explicitly correlated F12 methods

## Capabilities (CRITICAL)
- Ground-state electronic structure
- Geometry optimization and transition states
- Vibrational frequencies and thermochemistry
- Reaction pathways (NEB-TS, IRC)
- GOAT global optimizer and conformer generator
- Spectroscopic properties:
  - UV-Vis and fluorescence
  - IR and Raman (analytical intensities)
  - NMR chemical shifts and coupling constants
  - EPR/ESR (g-tensors, hyperfine coupling, ZFS)
  - Mössbauer parameters
  - VCD (vibrational circular dichroism)
  - X-ray absorption spectroscopy (XAS, XES)
- Magnetic properties at high-level ab initio methods
- Excited state calculations (TDDFT, EOM-CC)
- Multireference calculations for complex systems
- DLPNO methods for large molecules (100+ atoms)
- QM/MM calculations (QM/QM and QM/MM ONIOM)
- Periodic boundary conditions (limited)
- Automated composite methods
- Machine learning potential interfaces
- Explicit solvator module

**Sources**: Official ORCA documentation, cited in 7/7 source lists

## Key Strengths

### Spectroscopy Focus:
- Unmatched spectroscopic property calculations
- EPR, NMR, Mössbauer parameters
- UV-Vis with vibronic effects
- Resonance Raman spectra
- X-ray absorption and emission

### DLPNO Methods:
- Domain-Based Local Pair Natural Orbitals
- Near-linear scaling CCSD(T)
- Accurate for large molecules
- Benchmark quality at reduced cost

### User Experience:
- Simple keyword-based input
- Extensive documentation and tutorials
- Active user forum
- GaussView, Avogadro, Chemcraft compatibility
- Machine-readable property file output (v6.0+)

### Performance (v6.0+):
- Near-complete code rewrite
- Significant DFT speedups
- Improved NEB-TS for transition states
- Enhanced parallel efficiency

## Inputs & Outputs
- **Input formats**:
  - Simple input file (! keyword-based)
  - XYZ coordinate files
  - Flexible %block syntax
  - GBW checkpoint files
  
- **Output data types**:
  - Detailed output file (.out)
  - Energies, gradients, Hessians
  - Molecular orbitals (.gbw)
  - Spectroscopic data
  - Property file (.property.txt, JSON-like)

## Interfaces & Ecosystem
- **Visualization**:
  - orca_plot for orbitals and densities
  - Compatible with Avogadro, Chemcraft, UCSF ChimeraX
  
- **Utilities**:
  - orca_mapspc - spectrum processing
  - orca_vib - vibrational analysis
  - orca_casscf - CASSCF restart
  - Explicit solvator, micro-solvator
  
- **QM/MM**:
  - Interface to AMBER, CHARMM force fields
  - openCOSMO-RS interface
  
- **External Integrations**:
  - Machine learning potential interfaces
  - CFOUR, MRCC external interfaces
  - FMM for embedding


## Workflow and Usage

### Input Format:
ORCA uses a keyword-based free-format input style.

```
! PBE0 def2-TZVP Opt Freq
%pal nprocs 8 end

* xyz 0 1
  O  0.0 0.0 0.0
  H  0.0 0.7 0.0
  H  0.7 0.0 0.0
*
```

### Running ORCA:
```bash
# Requires full path usually
/path/to/orca input.inp > output.out
```

### Common Tasks:
- **Optimization**: `! Opt` keyword
- **Frequencies**: `! Freq` keyword
- **Excited States**: `%tdDFT` block
- **Spectroscopy**: `%eprnmr` block

## Advanced Features

### DLPNO-CCSD(T):
- Domain-Based Local Pair Natural Orbital method
- Near-linear scaling coupled cluster
- Enables "Gold Standard" accuracy for 100+ atoms
- Configurable accuracy settings (Loose, Normal, Tight)
- Major breakthrough for large systems

### Explicit Solvation:
- Automated placement of solvent molecules
- Micro-solvation environments
- Combined with continuum models (CPCM/SMD)
- Essential for spectroscopic accuracy

### Multireference Methods:
- CASSCF with automated starting orbitals
- NEVPT2 for dynamic correlation
- MRCI and MRCC support
- ROCIS for X-ray spectroscopy
- User-friendly "black box" approach

### GOAT Global Optimization:
- Global Optimization via Artificial Intelligence
- Conformer searching
- Reaction path finding
- Automated complex searches

### QM/MM Integration:
- ORCA/Gromacs interface
- Internal QM/MM module
- Electrostatic embedding
- Efficient for enzymatic capability

## Performance Characteristics
- **Speed**: RIJCOSX approximation speeds up Hybrid DFT by 10-100x
- **Scalability**: Good parallel scaling up to ~64-128 cores (varies by method)
- **Efficiency**: Optimized algorithms (integrals, SCF)
- **Memory**: Moderate requirements, disk usage can be high for coupled cluster
- **Disk I/O**: Heavy I/O for large calculations, fast SSDs recommended

## Computational Cost
- **DFT (RI)**: Very fast, linear scaling
- **Hybrid DFT**: Fast with RIJCOSX
- **MP2**: O(N^5), efficient with RI
- **DLPNO-CCSD(T)**: Near linear scaling, huge cost saving vs canonical
- **CASSCF**: Exponential scaling with active space
- **Canonical CCSD(T)**: O(N^7), very expensive

## Comparison with Other Codes
- **vs Gaussian**: ORCA free for academia, superior spectroscopy, comparable DFT speed; Gaussian has more composite methods.
- **vs NWChem**: ORCA easier input, better spectroscopy; NWChem scales better to 1000s of cores.
- **vs Turbomole**: Both very fast DFT; ORCA has broader multireference and spectroscopy features.
- **vs Dalton**: ORCA faster for standard properties; Dalton specialized for complex response theory.
- **Unique strength**: DLPNO methods for large-molecule accuracy, spectroscopy prediction, transition metal handling.

## Best Practices

### Basis Sets:
- Use "def2" family (def2-SVP, def2-TZVP, etc.)
- Always use auxiliary basis sets (`! def2/J` or `! def2/JK`) for RI speedup
- Check convergence for diffuse functions

### Approximation Techniques:
- Always use RI approximation for pure DFT
- Use RIJCOSX for hybrid functionals (standard in ORCA 5+)
- Use DLPNO for coupled cluster on >20 atoms

### Convergence:
- Use `! TightSCF` for spectroscopic properties
- Use `! SlowConv` for difficult transition metals
- Check spin contamination for open-shell systems

## Community and Support
- **Forum**: Very active ORCA Forum (crucial resource)
- **Manual**: Extremely detailed (>1000 pages) with theory reviews
- **Tutorials**: Library of official tutorials (Spring School)
- **Development**: Rapid release cycle (Max Planck Institute)
- **License**: Free for academic use, registration required

## Application Areas

### Transition Metal Chemistry:
- Open-shell systems
- Spin-state energetics
- Catalytic cycles
- Bioinorganic chemistry

### Spectroscopy:
- EPR/ESR of radicals and metal complexes
- NMR shieldings and coupling
- Mössbauer isomer shifts
- X-ray spectroscopy

### Large Molecules:
- DLPNO-CCSD(T) for 100+ atoms
- Protein-ligand binding energies
- Noncovalent interactions

## Limitations & Known Constraints
- **Registration required**: Free for academics but requires registration
- **Not open-source**: Binary distribution only
- **Periodic systems**: Limited support compared to solid-state codes
- **System size**: Molecular focus; not optimized for extended systems
- **Parallelization**: Excellent but varies by method
- **Platform**: Linux, macOS, Windows binaries
- **Documentation**: Extensive but requires quantum chemistry familiarity

## Verification & Sources
**Primary sources**:
1. Official website: https://orcaforum.kofo.mpg.de/
2. Manual: https://www.faccts.de/docs/orca/6.0/manual/
3. F. Neese, WIREs Comput. Mol. Sci. 12, e1606 (2022) - ORCA 5.0
4. F. Neese, WIREs Comput. Mol. Sci. 8, e1327 (2018) - Software update
5. ORCA 6.0 (July 2024), ORCA 6.1.0 (June 2025) - Latest releases

**Secondary sources**:
1. ORCA tutorials and workshops (Max Planck Institute)
2. Published applications in chemistry (>10,000 citations)
3. Spectroscopy benchmarks
4. Confirmed in 7/7 source lists (claude, g, gr, k, m, q, z)

**Confidence**: CONFIRMED - Appears in all 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE and ACCESSIBLE
- Software: Free academic binaries available
- Community support: Very active (forum, tutorials)
- Academic citations: >15,000 (various versions)
- Active development: Major v6.0 release (2024)
- Specialized strength: Spectroscopy, transition metals, DLPNO-CCSD(T), user-friendly
