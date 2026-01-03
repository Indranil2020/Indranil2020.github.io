# DP (Dielectric Properties)

## Official Resources
- Homepage: http://dp-code.org/
- Documentation: http://dp-code.org/documentation.html
- Source Repository: Distributed via website (GPL)
- License: GNU General Public License v2.0

## Overview
DP is a code for calculating the linear response dielectric properties of periodic systems (bulk crystals, surfaces, etc.) using Time-Dependent Density Functional Theory (TDDFT) in the linear response regime. It works in the frequency domain and uses a plane-wave basis set. DP calculates the macroscopic dielectric function, electron energy loss spectra (EELS), and IXS spectra.

**Scientific domain**: Dielectric properties, TDDFT, EELS, IXS  
**Target user community**: Spectroscopists, electronic structure theorists

## Theoretical Methods
- Time-Dependent Density Functional Theory (TDDFT)
- Linear Response Theory
- Random Phase Approximation (RPA)
- Adiabatic Local Density Approximation (ALDA)
- Long-range Component (LRC) kernels
- Plane-wave basis set implementation

## Capabilities (CRITICAL)
- Calculation of macroscopic dielectric function $\epsilon(\omega)$
- Electron Energy Loss Spectra (EELS) for low and high momentum transfer
- Inelastic X-ray Scattering (IXS) spectra
- Local field effects (LFE)
- Crystal local field effects
- Anisotropic dielectric tensors
- Integration with ABINIT (uses WFK files)

**Sources**: DP website, Comp. Phys. Comm. 144, 237 (2002)

## Inputs & Outputs
- **Input formats**: `dp.in` (input parameters), `WFK` (wavefunctions from ABINIT)
- **Output data types**: `out.eps` (dielectric function), `out.eels` (EELS spectrum)

## Interfaces & Ecosystem
- **ABINIT**: Primary interface (reads ABINIT binary wavefunction files)
- **Visualization**: Output is text-based for plotting (gnuplot, xmgrace)

## Workflow and Usage
1. Run ABINIT ground state calculation to generate wavefunctions (`_WFK`).
2. Create `dp.in`: Define frequency range, momentum transfer q, kernel.
3. Run `dp`: `dp < dp.in > dp.out`
4. Plot the resulting spectra.

## Performance Characteristics
- Frequency-domain calculation
- Scales with number of bands and plane waves
- Parallelization (MPI) supported

## Application Areas
- Optical absorption of solids
- Plasmon excitations
- Core-level spectroscopy (EELS)
- Dielectric screening

## Community and Support
- Developed by ETSF (European Theoretical Spectroscopy Facility) members
- Active user forum
- Long-standing code in the TDDFT community

## Verification & Sources
**Primary sources**:
1. Homepage: http://dp-code.org/
2. Publication: V. Olevano, L. Reining, and G. Siny, "DP: a code for dielectric properties", unpublished/website.

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: AVAILABLE
- Source: OPEN (GPL)
- Development: ACTIVE (ETSF)
- Applications: Dielectric function, TDDFT, EELS, IXS
