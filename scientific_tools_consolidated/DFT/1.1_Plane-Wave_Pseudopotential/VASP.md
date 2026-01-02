# VASP

## Official Resources
- Homepage: https://www.vasp.at/
- Documentation: https://www.vasp.at/wiki/
- Source Repository: Proprietary (source available to licensees)
- License: Commercial and academic licenses required

## Overview
VASP (Vienna Ab initio Simulation Package) is a leading commercial plane-wave DFT code for performing quantum mechanical calculations of atomic and electronic structures. Known for its robustness, comprehensive features, extensive validation, and exceptional performance, VASP is the most widely used DFT code in materials science and computational chemistry, with particularly strong implementations of hybrid functionals, many-body methods, and advanced algorithms.

**Scientific domain**: Plane-wave DFT, electronic structure, materials science, solid-state physics  
**Target user community**: Materials scientists, computational physicists, industrial researchers

## Theoretical Methods
- Kohn-Sham DFT (LDA, GGA, meta-GGA)
- Plane-wave basis with pseudopotentials
- Projector augmented wave (PAW) method
- Ultrasoft pseudopotentials (USPP)
- Hybrid functionals (PBE0, HSE06, B3LYP, etc.)
- Range-separated hybrids
- van der Waals corrections (DFT-D2/D3/D4, vdW-DF, TS, MBD)
- DFT+U for correlated systems
- Hartree-Fock
- GW approximation (G₀W₀, scGW, scGW₀)
- Random Phase Approximation (RPA)
- Bethe-Salpeter equation (BSE)
- Møller-Plesset perturbation theory (MP2)
- Coupled Cluster (CCSD, CCSD(T))
- Time-Dependent DFT (TDDFT)
- Density Functional Perturbation Theory (DFPT)
- Non-collinear magnetism
- Spin-orbit coupling
- Bethe-Salpeter Equation (BSE)
- Density Functional Perturbation Theory (DFPT) at Γ-point

## Capabilities (CRITICAL)
- Ground-state electronic structure calculations
- Geometry optimization (ionic relaxation)
- Transition state searches (NEB, dimer method)
- Ab initio molecular dynamics (AIMD): NVE, NVT, NPT ensembles
- Phonon calculations via finite differences or DFPT at Γ
- Elastic constants and mechanical properties
- Optical properties via frequency-dependent dielectric function
- Band structure and density of states
- Charge density and electrostatic potential analysis
- Magnetic properties (spin-polarized, non-collinear magnetism, spin-orbit coupling)
- Core-level spectroscopy (XPS, EELS)
- STM image simulation
- Constrained DFT calculations
- Machine learning force fields (VASP ML-FF)

**Sources**: Official VASP wiki, user manual (vasp.at/wiki), cited in 7/7 source lists

## Inputs & Outputs
- **Input formats**: 
  - POSCAR (structure file)
  - INCAR (control parameters)
  - POTCAR (pseudopotential/PAW data)
  - KPOINTS (k-point sampling)
  
- **Output data types**:
  - OUTCAR (main output, energies, forces, stresses)
  - CONTCAR (relaxed structure)
  - vasprun.xml (structured XML output)
  - CHGCAR, LOCPOT (charge densities, potentials)
  - EIGENVAL, DOSCAR (eigenvalues, DOS)
  - WAVECAR (wavefunctions, for restarts)

## Interfaces & Ecosystem
- **Framework integrations**:
  - ASE (Atomic Simulation Environment) - calculator interface
  - pymatgen - structure I/O, analysis
  - AiiDA - workflow automation (aiida-vasp plugin)
  - atomate/atomate2 - high-throughput workflows
  - custodian - error handling
  
- **Post-processing compatibility**:
  - vaspkit - comprehensive post-processing
  - sumo - band structure/DOS plotting
  - pyprocar - band structure analysis
  - Lobster - chemical bonding analysis
  - BoltzTraP - transport properties
  - Phonopy - phonon calculations
  - Wannier90 - Wannier function generation
  
- **Machine learning interfaces**:
  - On-the-fly ML potential training (VASP 6.x)
  - Integration with LAMMPS via ML-FF

## Limitations & Known Constraints
- **Licensing**: Requires commercial or academic license; not open-source
- **Cost**: Significant license fees for commercial use
- **Scaling**: Good parallel scaling up to ~1000 cores for large systems; efficiency depends on system size and k-point sampling
- **System size**: Practical limit ~500-1000 atoms for standard DFT; larger systems possible with ML-FF
- **Pseudopotentials**: PAW datasets must be obtained separately; quality depends on PAW library version
- **Convergence**: Can be challenging for metallic systems or strongly correlated materials without appropriate settings
- **Hybrid functionals**: Computationally expensive; limited to smaller systems (~100 atoms)
- **Documentation**: Some advanced features poorly documented; community wiki varies in quality

## Verification & Sources
**Primary sources**:
1. Official VASP website: https://www.vasp.at/
2. VASP wiki documentation: https://www.vasp.at/wiki/
3. G. Kresse and J. Furthmüller, Phys. Rev. B 54, 11169 (1996) - VASP original paper
4. G. Kresse and D. Joubert, Phys. Rev. B 59, 1758 (1999) - PAW implementation

**Secondary sources**:
1. ASE calculator documentation: https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html
2. pymatgen VASP I/O: https://pymatgen.org/pymatgen.io.vasp.html
3. Materials Project workflows (atomate): Uses VASP as primary DFT engine
4. Confirmed in 7/7 source lists (claude, g, gr, k, m, q, z)

**Confidence**: CONFIRMED - Appears in all 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: ACCESSIBLE
- Community support: EXTENSIVE (VASP forum, Materials Project, stackoverflow)
- Academic citations: >40,000 (Google Scholar)
