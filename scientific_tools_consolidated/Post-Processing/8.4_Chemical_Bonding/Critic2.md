# Critic2

## Official Resources
- Homepage: https://github.com/aoterodelaroza/critic2
- Documentation: https://aoterodelaroza.github.io/critic2/
- Source Repository: https://github.com/aoterodelaroza/critic2
- License: GNU General Public License v3.0

## Overview
Critic2 is a code for the analysis of quantum mechanical electron density and other scalar fields in periodic solids. It implements the Quantum Theory of Atoms in Molecules (QTAIM) and other topological analysis methods (NCI, ELF, etc.) for crystals. It can read output from a wide variety of electronic structure codes (WIEN2k, VASP, QE, Abinit, etc.) and perform critical point finding, basin integration, and chemical bonding analysis.

**Scientific domain**: QTAIM, electron density topology, chemical bonding  
**Target user community**: Solid-state chemists, crystallographers, DFT users

## Theoretical Methods
- Quantum Theory of Atoms in Molecules (QTAIM)
- Critical point (CP) search (Morse theory)
- Basin integration (atomic charges/volumes)
- Non-Covalent Interactions (NCI) index
- Electron Localization Function (ELF)
- Laplacian of electron density
- Interacting Quantum Atoms (IQA) - implementation dependent

## Capabilities (CRITICAL)
- Finding critical points (bond, ring, cage, nuclear) in electron density
- Integration of properties over atomic basins (charges, volumes, dipole moments)
- NCI plot analysis for weak interactions
- ELF/ELI analysis
- Analysis of VASP, WIEN2k, Quantum ESPRESSO, Abinit, Elk, CP2K, Gaussian (cube) output
- Handling of periodic boundary conditions
- Automatic chemical bond detection

**Sources**: Critic2 documentation, Comp. Phys. Comm. 185, 1007 (2014)

## Key Strengths

### Multi-Code Support:
- WIEN2k, VASP, QE, Abinit
- Gaussian, ORCA cube files
- Elk, CP2K support
- Consistent analysis

### Comprehensive Analysis:
- QTAIM critical points
- NCI visualization
- ELF/ELI analysis
- Basin integration

### Open Source:
- GPL v3 licensed
- Active development
- GitHub hosted
- Extensive documentation

## Inputs & Outputs
- **Input formats**: Code-specific density files (e.g., CHGCAR, rho.r, cube), Critic2 input script (.cri)
- **Output data types**: Critical point list, integrated basin properties, visualization files (cube, xyz, vtk)

## Interfaces & Ecosystem
- **Supported Codes**: WIEN2k, VASP, QE, Abinit, Elk, CP2K, Gaussian, Orca, etc.
- **Visualization**: Generates output for Jmol, VMD, Vest, gnuplot
- **Scripting**: Command-line interface with scripting capabilities

## Workflow and Usage
1. Perform DFT calculation and write electron density.
2. Prepare Critic2 input: Define crystal structure and density file.
3. Run Critic2: `critic2 input.cri`
4. Common tasks: `auto` (find all CPs), `basin` (integrate basins), `nci` (compute NCI).

## Performance Characteristics
- Efficient critical point search
- Basin integration can be computationally intensive (grid-based)
- OpenMP parallelization

## Limitations & Known Constraints
- **Learning curve**: Requires QTAIM knowledge
- **Basin integration**: Can be computationally intensive
- **Grid resolution**: Trade-off with accuracy
- **Input preparation**: Code-specific density formats

## Comparison with Other Tools
- **vs AIMAll**: Critic2 open-source, AIMAll commercial
- **vs Bader (Henkelman)**: Critic2 more features, Bader faster
- **vs Multiwfn**: Critic2 periodic systems, Multiwfn molecules
- **Unique strength**: Multi-code support, periodic QTAIM

## Application Areas
- Bonding characterization in solids
- Hydrogen bonding and weak interactions (NCI)
- High-pressure phases
- Charge transfer analysis
- Pore analysis in MOFs/Zeolites

## Best Practices
- Use converged DFT density
- Choose appropriate grid for basin integration
- Validate critical points with chemical intuition
- Cross-check with NCIPLOT for weak interactions

## Community and Support
- Open-source (GPL v3)
- Developed by A. Otero-de-la-Roza (University of Oviedo)
- Active GitHub repository

## Verification & Sources
**Primary sources**:
1. Homepage: https://github.com/aoterodelaroza/critic2
2. Documentation: https://aoterodelaroza.github.io/critic2/
3. Publication: A. Otero-de-la-Roza et al., Comp. Phys. Comm. 185, 1007 (2014)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN (GitHub)
- Development: ACTIVE (Otero-de-la-Roza)
- Applications: QTAIM, NCI, critical points, solids
