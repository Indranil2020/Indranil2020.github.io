# GAPW (Gaussian and Augmented Plane Waves)

## **METHOD NOT SOFTWARE**

## Clarification
**GAPW** is **NOT** a standalone software package. It is a **method** (Gaussian and Augmented Plane Waves) used for electronic structure calculations.

## Description
The GAPW method uses a dual basis set approach:
- **Gaussian functions** for describing the atomic orbitals and the core region.
- **Augmented Plane Waves** (or auxiliary plane waves) for the smooth part of the density and potential.

This method combines the efficiency of plane waves with the compactness of Gaussian functions, allowing for all-electron calculations or accurate pseudopotential calculations.

## Implementation
The GAPW method is most notably implemented in:
- **CP2K**: The Quickstep module of CP2K is the primary reference implementation for GAPW.
- **GPAW**: While primarily known for PAW, methods combining localized and delocalized functions share theoretical grounds.

## Recommendation
If you are looking for software that implements the GAPW method, please refer to:
- **CP2K** (Category 1.1)

## Verification
- **Primary Source**: CP2K Documentation (https://www.cp2k.org/)
- **Literature**: Lippert, G., Hutter, J., & Parrinello, M. (1997). "A hybrid Gaussian and plane wave density functional scheme". Molecular Physics.
