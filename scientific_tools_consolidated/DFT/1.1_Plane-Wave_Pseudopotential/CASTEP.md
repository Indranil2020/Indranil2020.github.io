# CASTEP

## Official Resources
- Homepage: http://www.castep.org/
- Documentation: http://www.castep.org/CASTEP/Documentation
- Source Repository: Proprietary (source available to licensees)
- License: Academic and commercial licenses available

## Overview
CASTEP is a leading academic and commercial plane-wave DFT code for studying materials from first principles. Developed in the UK, it provides comprehensive capabilities for calculating properties of materials including metals, semiconductors, ceramics, and molecular systems, with particular strengths in spectroscopy calculations (NMR, EPR, optical) and density functional perturbation theory.

**Scientific domain**: Plane-wave DFT, materials science, spectroscopy, solid-state physics  
**Target user community**: Materials scientists, solid-state physicists, spectroscopists

## Theoretical Methods
- Kohn-Sham DFT (LDA, GGA, meta-GGA)
- Hybrid functionals (PBE0, HSE06, B3LYP)
- Plane-wave basis with pseudopotentials
- Ultrasoft pseudopotentials (USPP)
- Projector augmented wave (PAW)
- Norm-conserving pseudopotentials
- Density Functional Perturbation Theory (DFPT)
- Linear response for phonons
- Electric field perturbations
- Magnetic field perturbations
- van der Waals corrections (DFT-D, TS, MBD)
- On-the-fly pseudopotential generation
- LDA, GGA, meta-GGA functionals
- Hybrid functionals (B3LYP, PBE0, HSE)
- DFT+U for correlated systems
- van der Waals corrections (DFT-D, TS)
- Spin-orbit coupling
- Non-collinear magnetism

## Capabilities (CRITICAL)
- Ground-state electronic structure
- Geometry optimization (BFGS, LBFGS, damped MD)
- Transition state searches (QST, NEB, dimer)
- Molecular dynamics (NVE, NVT, NPT, NPH)
- Phonon calculations via finite displacement or DFPT
- Elastic and mechanical properties
- Dielectric properties and Born effective charges
- Optical properties (absorption, reflectivity)
- NMR chemical shifts and J-coupling
- EPR g-tensors and hyperfine tensors
- Raman and IR intensities
- Core-level spectroscopy (XPS, EELS)
- Electric field gradients
- Band structure and density of states
- Wannier functions
- STM/AFM image simulation
- Accurate stress tensor calculations
- Equation of state fitting

**Sources**: Official CASTEP documentation, cited in 6/7 source lists

## Inputs & Outputs
- **Input formats**:
  - .cell file (unit cell and atomic positions)
  - .param file (calculation parameters)
  - .check file (checkpoints for restarts)
  - Standard structure formats via conversion
  
- **Output data types**:
  - .castep file (main output with energies, forces)
  - .geom file (optimized geometries)
  - .md file (molecular dynamics trajectories)
  - .phonon file (phonon frequencies and eigenvectors)
  - .bands file (electronic band structures)
  - .pdos file (projected density of states)
  - .cst_esp file (electrostatic potential)

## Interfaces & Ecosystem
- **Framework integrations**:
  - ASE - calculator interface
  - pymatgen - structure I/O and analysis
  - Materials Studio - GUI interface (commercial)
  - Phonopy - phonon post-processing
  
- **Pre/Post-processing**:
  - c2x - CASTEP to various formats converter
  - OptaDOS - advanced DOS analysis
  - dos.pl - DOS plotting utility
  - dispersion.pl - phonon dispersion plotting
  
- **Workflow integration**:
  - Can be integrated into high-throughput workflows
  - Grid computing support

## Limitations & Known Constraints
- **Licensing**: Requires academic or commercial license; not open-source
- **Cost**: License fees for commercial use
- **Pseudopotentials**: Quality depends on pseudopotential library used
- **Memory**: Plane-wave methods memory-intensive for large systems
- **Hybrid functionals**: Computationally expensive; limited to smaller systems
- **System size**: Practical limit ~500-1000 atoms for standard DFT
- **Convergence**: Metals require careful smearing parameter selection
- **Documentation**: Comprehensive but requires familiarity with input file format
- **Platform support**: Primarily Linux/Unix; Windows support limited

## Verification & Sources
**Primary sources**:
1. Official website: https://www.castep.org/
2. Documentation: http://www.castep.org/CASTEP/Documentation
3. S. J. Clark et al., Z. Kristallogr. 220, 567 (2005) - CASTEP first principles code
4. M. D. Segall et al., J. Phys.: Condens. Matter 14, 2717 (2002) - CASTEP method

**Secondary sources**:
1. CASTEP tutorials and workshops
2. ASE calculator documentation
3. Materials Studio documentation (GUI for CASTEP)
4. Confirmed in 6/7 source lists (claude, g, gr, k, m, q)

**Confidence**: CONFIRMED - Appears in 6 of 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: ACCESSIBLE
- License: Academic/Commercial (verified)
- Community support: Active (user forums, support)
- Academic citations: >3,000 (main papers)
- Industrial use: Extensive in materials science
