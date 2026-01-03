# FEFF

## Official Resources
- Homepage: https://feff.uw.edu/
- Documentation: https://feff.uw.edu/documentation/
- Source Repository: Proprietary (Academic/Commercial licenses)
- License: Proprietary

## Overview
FEFF is an automated program for ab initio multiple scattering calculations of X-ray Absorption Fine Structure (XAFS), X-ray Absorption Near-Edge Structure (XANES), and various other spectroscopies for clusters of atoms. Developed at the University of Washington, FEFF uses a real-space Green's function approach, making it highly effective for non-periodic systems, nanoparticles, and defects, as well as crystals.

**Scientific domain**: X-ray spectroscopy (XAFS, XANES, EXAFS), multiple scattering theory  
**Target user community**: Spectroscopists, materials scientists, chemists

## Theoretical Methods
- Real-space Multiple Scattering Theory (RSMS)
- Green's function formalism
- Self-consistent field (SCF) potentials
- Muff-tin approximation
- Full multiple scattering (FMS)
- Time-dependent DFT (TDDFT) for core-hole screening (in FEFF9)
- Many-pole self-energy (GW approximation)

## Capabilities (CRITICAL)
- Calculation of EXAFS and XANES spectra
- Electron Energy Loss Spectroscopy (EELS)
- X-ray Emission Spectroscopy (XES)
- Compton scattering
- Local Density of States (LDOS)
- X-ray Magnetic Circular Dichroism (XMCD)
- Non-resonant Inelastic X-ray Scattering (NRIXS)
- Core-hole effects
- Debye-Waller factors via correlated Debye model

**Sources**: FEFF website, Rev. Mod. Phys. 72, 621 (2000)

## Inputs & Outputs
- **Input formats**: `feff.inp` (geometry, potentials, control flags)
- **Output data types**: `xmu.dat` (absorption cross-section), `chi.dat` (EXAFS), `ldos.dat`, `feff.bin`

## Interfaces & Ecosystem
- **JFEFF**: GUI for input generation
- **Athena/Artemis (Demeter)**: Standard analysis software that uses FEFF for fitting
- **Larch**: Analysis library wrapping FEFF
- **ASE**: Interface available

## Workflow and Usage
1. Create `feff.inp` file (atomic coordinates, potentials).
2. Run FEFF modules:
   - `pot`: Calculate potentials
   - `xsph`: Phase shifts
   - `fms`: Full multiple scattering
   - `path`: Path expansion (for EXAFS)
   - `genfmt`: XAFS parameters
   - `ff2chi`: Chi calculation
3. Analyze `xmu.dat` or use output for fitting in Artemis.

## Performance Characteristics
- Highly efficient for high-energy EXAFS (path expansion)
- XANES (FMS) scales with cluster size (up to hundreds of atoms)
- Parallelized (MPI)

## Application Areas
- Structure determination from EXAFS (coordination numbers, bond lengths)
- Catalyst characterization (nanoparticles)
- Biological metalloproteins
- Amorphous materials and liquids
- Phase identification via XANES

## Community and Support
- Developed by Rehr Group (University of Washington)
- Very large user base in spectroscopy community
- Commercial support available
- Extensive workshops and training

## Verification & Sources
**Primary sources**:
1. Homepage: https://feff.uw.edu/
2. Publication: J. J. Rehr and R. C. Albers, Rev. Mod. Phys. 72, 621 (2000)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: PROPRIETARY (Widely used)
- Development: ACTIVE (Rehr Group)
- Applications: XAFS, XANES, EELS, multiple scattering
