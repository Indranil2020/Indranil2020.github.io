# FDMNES

## Official Resources
- Homepage: http://neel.cnrs.fr/fdmnes
- Documentation: http://neel.cnrs.fr/fdmnes/Manual/Manual.html
- Source Repository: Distributed via website (Open source)
- License: GNU General Public License v3.0 (since 2020)

## Overview
FDMNES is a code for simulating X-ray Absorption Spectroscopy (XAS) including XANES and EXAFS, X-ray Emission Spectroscopy (XES), and Resonant Inelastic X-ray Scattering (RIXS). It employs both Finite Difference Method (FDM) and Multiple Scattering Theory (MST) to calculate excited states in clusters or periodic systems. It is particularly powerful for handling low-symmetry systems and core-hole effects.

**Scientific domain**: X-ray spectroscopy, multiple scattering, finite difference method  
**Target user community**: Spectroscopists, materials scientists, beamline scientists

## Theoretical Methods
- Finite Difference Method (FDM) solving Schrödinger equation
- Multiple Scattering Theory (MST) / Green's function (Muff-tin)
- Time-Dependent DFT (TDDFT) for core-hole screening
- Spin-orbit coupling (fully relativistic)
- Hubbard U corrections (LSDA+U)
- Tensor calculations for dichroism (XMCD, XNCD)

## Capabilities (CRITICAL)
- Calculation of K, L, M edges XANES and EXAFS
- RIXS and XES simulations
- Choice between FDM (more accurate, slower) and MST (faster, muffin-tin approx)
- Non-muffin-tin potentials (FDM)
- Spin-polarized calculations
- Analysis of dichroism and optical activity
- GUI for input generation and visualization

**Sources**: FDMNES website, Phys. Rev. B 63, 125120 (2001)

## Key Strengths

### Dual Methods:
- FDM for accuracy
- MST for speed
- User choice
- Non-muffin-tin option

### Comprehensive Features:
- XANES, EXAFS, RIXS
- Dichroism (XMCD, XNCD)
- Spin-orbit coupling
- GUI available

### Active Development:
- GPL licensed
- Regular updates
- Workshops available
- Good documentation

## Inputs & Outputs
- **Input formats**: `fdmnes.ind` (input parameters), `fdmfile.txt` (structure)
- **Output data types**: `xanes.txt` (spectra), `density.txt` (DOS), `bse.txt`

## Interfaces & Ecosystem
- **GUI**: Java-based graphical interface provided
- **Fit**: Integration with fitting procedures
- **Parallelization**: MPI support for large calculations

## Workflow and Usage
1. Define structure (atoms, space group).
2. Choose method (Green/MST or FDM).
3. Set parameters (radius, edges, polarization).
4. Run `fdmnes`.
5. Convolution: Broaden raw spectra with core-hole lifetime.

## Performance Characteristics
- MST is very fast (seconds/minutes)
- FDM is computationally intensive but handles non-spherical potentials accurately
- Memory scaling with FDM grid size

## Limitations & Known Constraints
- **FDM cost**: Slower than MST
- **Learning curve**: Many parameters
- **Grid convergence**: FDM needs testing
- **Documentation**: Could be more extensive

## Comparison with Other Tools
- **vs FEFF**: FDMNES FDM option, FEFF muffin-tin only
- **vs xspectra**: Different theoretical approaches
- **vs OCEAN**: FDMNES faster, OCEAN BSE-based
- **Unique strength**: FDM for non-spherical potentials

## Application Areas
- Structure refinement from XANES
- Magnetic materials (XMCD)
- Biological metalloproteins
- Distorted local environments
- Actinides and lanthanides (relativistic effects)

## Best Practices
- Start with MST for quick tests
- Use FDM for final accuracy
- Converge cluster size
- Apply appropriate broadening

## Community and Support
- Developed at Institut Néel (CNRS/Grenoble)
- Active mailing list
- Frequent updates and workshops

## Verification & Sources
**Primary sources**:
1. Homepage: http://neel.cnrs.fr/fdmnes
2. Publication: Y. Joly, Phys. Rev. B 63, 125120 (2001)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: AVAILABLE
- Source: OPEN (GPL v3)
- Development: ACTIVE (Y. Joly, O. Bunau)
- Applications: XAS, FDM, MST, spectroscopy
