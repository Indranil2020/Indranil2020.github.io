# Quantum ESPRESSO

## Official Resources
- Homepage: https://www.quantum-espresso.org/
- Documentation: https://www.quantum-espresso.org/documentation/
- Source Repository: https://gitlab.com/QEF/q-e
- License: GNU General Public License v2.0

## Overview
Quantum ESPRESSO (QE) is an integrated suite of open-source codes for electronic-structure calculations and materials modeling at the nanoscale based on density-functional theory, plane waves, and pseudopotentials. It is one of the most widely used DFT codes worldwide, known for its comprehensive features, excellent parallelization, and strong community support, particularly for solid-state and materials science applications.

**Scientific domain**: Plane-wave DFT, electronic structure, materials science, solid-state physics  
**Target user community**: Materials scientists, condensed matter physicists, computational chemists

## Theoretical Methods
- Kohn-Sham DFT (LDA, GGA, meta-GGA)
- Plane-wave basis with pseudopotentials
- Ultrasoft pseudopotentials (USPP)
- Projector augmented wave (PAW)
- Norm-conserving pseudopotentials
- Density Functional Perturbation Theory (DFPT)
- Time-Dependent DFT (TDDFT, turbo_lanczos)
- Hybrid functionals (via ACE, SternheimerGW)
- van der Waals corrections (vdW-DF, DFT-D2/D3, Tkatchenko-Scheffler)
- DFT+U for correlated systems
- Non-collinear magnetism
- Spin-orbit coupling
- GW approximation (via SternheimerGW)
- Quantum transport (via Wannier90+WanT)

## Capabilities (CRITICAL)
- Ground-state electronic structure (pw.x module)
- Structure optimization and variable-cell relaxation
- Phonon calculations via DFPT (ph.x module)
- Electron-phonon coupling (EPW code)
- Ab initio molecular dynamics: Born-Oppenheimer and Car-Parrinello
- Linear-response TDDFT (turbo_lanczos.x, turbo_davidson.x)
- X-ray absorption spectroscopy (xspectra.x)
- NMR and EPR parameters
- STM image simulation
- Wannier function generation (interface to Wannier90)
- Transport properties via NEGF (pwcond.x)
- GW and BSE calculations (via WEST or yambo post-processing)
- Post-processing: bands, DOS, charge densities, potentials
- Constrained DFT
- Hubbard corrections (DFT+U)

**Sources**: Official QE website, user guide, cited in 7/7 source lists

## Inputs & Outputs
- **Input formats**:
  - Main input file (pw.in) with Fortran namelist format
  - Pseudopotential files (.UPF format)
  - Atomic positions (crystal, alat, bohr, angstrom)
  - K-point specifications
  
- **Output data types**:
  - Standard output (text format with energies, forces, convergence)
  - XML data file (data-file-schema.xml)
  - Charge density (charge-density.dat)
  - Wavefunction files (wfc*.dat)
  - Phonon dynamical matrices
  - Force constants for phonons
  - XML outputs for post-processing

## Interfaces & Ecosystem
- **Framework integrations**:
  - ASE - native calculator support
  - AiiDA - comprehensive plugin (aiida-quantumespresso)
  - pymatgen - structure I/O and analysis
  - atomate2 - workflow support
  - Phonopy - phonon post-processing
  - Wannier90 - tight integration via pw2wannier90.x
  
- **Post-processing tools**:
  - Integrated: pp.x, bands.x, dos.x, projwfc.x, plotband.x
  - External: yambo (GW/BSE), BerkeleyGW, EPW (e-ph)
  - Python interfaces: qe-tools, aiida-quantumespresso
  
- **Workflow integration**:
  - AiiDA workflows for phonons, bands, convergence studies
  - Materials Project workflows via atomate

## Limitations & Known Constraints
- **Scaling**: Good parallel scaling up to thousands of cores for large systems; k-point parallelization most efficient
- **System size**: Practical limit ~500-1000 atoms for standard plane-wave DFT; larger with specific optimizations
- **Pseudopotentials**: Quality depends on pseudopotential library; requires careful validation
- **Memory**: Can be memory-intensive for hybrid functionals or exact exchange
- **Convergence**: Metallic systems require appropriate smearing; insulators need careful k-point sampling
- **Hybrid functionals**: Computationally expensive via experimental ACE feature

## Computational Cost
- **Plane Wave DFT**: $O(N^3)$ with system size.
- **Parallel Scaling**: Excellent scaling (up to thousands of cores) due to FFT parallelization.
- **Memory**: Plane-wave basis requires significant memory for vacuum regions (use implicit solvent or Coulomb cutoff if needed).
- **User interface**: Command-line driven; learning curve for input file syntax
- **Documentation**: Extensive but can be overwhelming; scattered across multiple manuals

## Comparison with Other Codes
- **vs VASP**: QE is free/open-source; performance is competitive (often slightly slower per SCF, but scales better to massive core counts).
- **vs ABINIT**: QE is generally faster for ground state/MD; ABINIT has historically stronger many-body perturbation theory features (though QE's Yambo interface is excellent).
- **vs CP2K**: QE uses pure plane waves (better for accuracy/hard potentials); CP2K uses Gaussian+Plane Waves (faster for large systems/MD).

## Best Practices
- **Parallelization**: Use hybrid MPI/OpenMP. Segregate cores with `-ndiag` for large matrices.
- **Smearing**: Use `smearing='mv'` (Marzari-Vanderbilt) for metals.
- **Pseudopotentials**: Stick to `SSSP` (Standard Solid State Pseudopotentials) library for verified accuracy.
- **Convergence**: Checking `ecutwfc` and `ecutrho` is critical; default values often insufficient for specialized properties.

## Community and Support
- **Forum**: Active discourse forum and mailing list.
- **Events**: Regular schools in Trieste (ICTP) and worldwide.
- **Development**: Open contributions via GitLab.

## Verification & Sources
**Primary sources**:
1. Official website: https://www.quantum-espresso.org/
2. User guide: https://www.quantum-espresso.org/Doc/user_guide/
3. P. Giannozzi et al., J. Phys.: Condens. Matter 21, 395502 (2009) - QE original paper
4. P. Giannozzi et al., J. Phys.: Condens. Matter 29, 465901 (2017) - QE update
5. P. Giannozzi et al., J. Chem. Phys. 152, 154105 (2020) - QE 6.5 features

**Secondary sources**:
1. ASE calculator: https://wiki.fysik.dtu.dk/ase/ase/calculators/espresso.html
2. AiiDA-QuantumESPRESSO plugin: https://github.com/aiidateam/aiida-quantumespresso
3. EPW documentation: https://docs.epw-code.org/
4. Confirmed in 7/7 source lists (claude, g, gr, k, m, q, z)

**Confidence**: CONFIRMED - Appears in all 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: ACCESSIBLE
- Source code: OPEN (GitLab)
- Community support: EXTENSIVE (user forum, mailing list)
- Academic citations: >15,000 (Google Scholar)
