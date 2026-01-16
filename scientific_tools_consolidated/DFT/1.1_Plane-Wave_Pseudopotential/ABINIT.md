# ABINIT

## Official Resources
- Homepage: https://www.abinit.org/
- Documentation: https://docs.abinit.org/
- Source Repository: https://github.com/abinit/abinit
- License: GNU General Public License v3.0

## Overview
ABINIT is a comprehensive open-source package for electronic structure calculations based on density-functional theory, using pseudopotentials and a plane-wave or wavelet basis set. It is particularly strong in linear-response calculations, many-body perturbation theory (GW), and excited-state methods.

**Scientific domain**: Condensed matter physics, materials science, electronic structure  
**Target user community**: Academic researchers, materials scientists requiring advanced response properties

## Theoretical Methods
- Density Functional Theory (DFT)
- Plane-wave pseudopotentials and PAW
- Wavelet basis sets (BigDFT integration)
- Density Functional Perturbation Theory (DFPT)
- Many-Body Perturbation Theory (GW approximation)
- Bethe-Salpeter Equation (BSE)
- Time-Dependent DFT (TDDFT)
- Dynamical Mean-Field Theory (DMFT)
- Constrained DFT
- Hybrid functionals
- DFT+U for correlated systems

## Capabilities (CRITICAL)
- Ground-state electronic structure calculations
- Geometry optimization and cell relaxation
- Molecular dynamics (Born-Oppenheimer, Langevin, Nosé-Hoover)
- Phonon calculations via DFPT (full q-point grid)
- Dielectric response and Born effective charges
- Elastic constants via response functions
- Electron-phonon coupling via DFPT
- GW quasiparticle energies (G₀W₀, self-consistent GW)
- BSE for optical absorption including excitonic effects
- TDDFT for excited states
- Non-linear optical properties
- Temperature-dependent electronic structure
- Spin-orbit coupling and non-collinear magnetism
- Electric field responses (Berry phase)
- Wannier function generation (interface to Wannier90)
- DFT+DMFT for strongly correlated systems
- PAW datasets generation

**Sources**: Official ABINIT documentation, tutorials, cited in 7/7 source lists

## Inputs & Outputs
- **Input formats**:
  - Main input file (.in or .abi) with structured keywords
  - Pseudopotential files (.psp8, .xml for PAW)
  - Density/wavefunction files for restarts
  - Structure files (various formats via ASE/pymatgen conversion)
  
- **Output data types**:
  - Main output (.out or .abo) with energies, forces, stresses
  - NetCDF files for density, wavefunctions, response functions
  - Phonon data (dynamical matrices, IFCs)
  - GW and BSE outputs
  - Formatted text files for band structures, DOS

## Interfaces & Ecosystem
- **Framework integrations**:
  - ASE - calculator interface
  - AiiDA - aiida-abinit plugin for workflows
  - pymatgen - structure I/O and analysis
  - Phonopy - phonon post-processing from force constants
  - Wannier90 - tight-binding Hamiltonian generation
  
- **Post-processing tools**:
  - abipy - Python package for ABINIT workflows and analysis
  - AbiPy post-processing utilities
  - cut3d - for visualizing densities and potentials
  - Built-in analysis tools (abicheck, abiview, etc.)
  
- **Related codes**:
  - BigDFT - wavelet basis integration
  - Yambo - alternative GW/BSE post-processor
  - BerkeleyGW - can read ABINIT wavefunctions

## Limitations & Known Constraints
- **Learning curve**: Input file syntax complex with many keywords; steep learning curve
- **Documentation**: Extensive but can be difficult to navigate; variable quality across topics
- **Performance**: Parallelization efficiency depends on calculation type; k-point parallelization most efficient
- **Memory**: GW and BSE calculations memory-intensive for large systems
- **Pseudopotentials**: Quality depends on pseudopotential table; requires validation
- **Convergence**: Response function calculations require careful convergence testing
- **Hybrid functionals**: Expensive; limited to smaller systems
- **Installation**: Build process can be complex with many optional dependencies
- **File I/O**: NetCDF files can become very large for big systems

## Computational Cost
- **Ground State (PW)**: $O(N^3)$.
- **GW/BSE**: Highly expensive ($O(N^4)$); requires massive memory.
- **Wavelets (BigDFT)**: $O(N)$ linear scaling mode available via BigDFT integration.

## Comparison with Other Codes
- **vs Quantum ESPRESSO**: ABINIT has a longer history with GW/response functions; QE is often faster for standard MD.
- **vs Yambo**: Yambo acts as a post-processor for QE; ABINIT has GW built-in, offering a more unified but sometimes monolithic experience.
- **vs VASP**: ABINIT is open-source and has more diverse basis set options (wavelets); VASP is faster for simple relaxation.

## Best Practices
- **Parallelization**: Study the `autoparal` feature (automatic parallelization tuning).
- **Pseudopotentials**: Use `PseudoDojo` tables (standard for ABINIT).
- **Convergence**: GW calculations require convergence of empty states (`nband`), which is much harder than ground state.

## Community and Support
- **Forum**: Official ABINIT Forum (forum.abinit.org).
- **Events**: ABINIT schools held annually (often in Europe/Louvain).
- **Development**: Hosted on GitHub (switched from diverse repos).

## Verification & Sources
**Primary sources**:
1. Official website: https://www.abinit.org/
2. Documentation: https://docs.abinit.org/
3. X. Gonze et al., Comput. Phys. Commun. 180, 2582 (2009) - ABINIT overview
4. X. Gonze et al., Comput. Mater. Sci. 25, 478 (2002) - ABINIT first paper
5. F. Bottin et al., Comput. Mater. Sci. 42, 329 (2008) - PAW implementation

**Secondary sources**:
1. ABINIT tutorials: https://docs.abinit.org/tutorial/
2. abipy documentation: http://abinit.github.io/abipy/
3. ASE calculator: https://wiki.fysik.dtu.dk/ase/ase/calculators/abinit.html
4. Confirmed in 7/7 source lists (claude, g, gr, k, m, q, z)

**Confidence**: CONFIRMED - Appears in all 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE and ACCESSIBLE
- Source code: OPEN (GitHub)
- Community support: Active (forum, mailing list)
- Academic citations: >3,000 (main papers)
- Active development: Regular releases
