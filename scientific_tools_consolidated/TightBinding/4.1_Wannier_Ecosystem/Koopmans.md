# Koopmans

## Official Resources
- **Homepage**: https://koopmans-functionals.org/
- **Documentation**: https://koopmans.readthedocs.io/
- **Source Repository**: https://github.com/koopmans-functionals/koopmans
- **License**: GNU General Public License v3.0

## Overview
**Koopmans** is a spectral functional code designed to predict accurate electronic structures, particularly band gaps and ionization potentials, by enforcing the **Generalized Koopmans' Theorem (GKT)**. It serves as a wrapper and extension around **Quantum ESPRESSO**, using **Wannier90** to construct localized orbitals on which Koopmans-compliant corrections are applied. The code addresses the band gap problem in standard DFT by ensuring that orbital energies correspond directly to charged excitation energies, achieving accuracy comparable to high-level GW calculations but at a significantly lower computational cost.

**Scientific domain**: Band structure theory, Spectroscopy, Insulators and Semiconductors, Photovoltaics
**Target user community**: Researchers needing accurate band gaps and spectral properties without the cost of GW

## Theoretical Methods
- **Generalized Koopmans' Theorem (GKT)**: Enforces linear behavior of energy with respect to fractional particle number.
- **Koopmans-Compliant Functionals**:
  - **KI (Koopmans Integral)**: Removes curvature in energy vs. occupation.
  - **KIPZ (Koopmans Integral + Perdew-Zunger)**: Includes screened self-interaction correction.
  - **KIER (Koopmans Integral + External Relaxed)**: accounts for orbital relaxation.
- **Screening Parameters**: Calculated ab initio via linear response or predicted via Machine Learning.
- **Orbital Minimization**: Variational optimization of orbital densities.
- **Wannier Localization**: Uses Wannier functions (from Wannier90) as the variational basis for corrections.

## Capabilities
- **Spectral Properties**:
  - Accurate Band Gaps
  - Ionization Potentials
  - Electron Affinities
  - Photoemission Spectra
  - Band-edge alignments
- **Automated Workflows**: 
  - Automatic determination of screening parameters.
  - End-to-end calculation from DFT ground state to spectral properties.
- **Machine Learning Integration**: Accelerates screening parameter prediction for large systems.
- **System Types**:
  - Molecules (finite systems)
  - Crystals (periodic systems)
  - Disordered solids and liquids

## Key Strengths
- **Accuracy vs. Cost**: Delivers GW-level accuracy for band gaps at a computational cost comparable to standard DFT (or hybrid functionals).
- **Physical Insight**: Provides a clear orbital-based picture of excitations using Wannier functions.
- **Automation**: Streamlines the complex process of calculating screening and corrections.
- **Integration**: Seamlessly works with the widely used Quantum ESPRESSO and Wannier90 ecosystems.

## Inputs & Outputs
- **Inputs**:
  - Quantum ESPRESSO input files (pw.x)
  - Wannier90 input files (wannier90.win)
  - JSON-based workflow configuration
- **Outputs**:
  - Corrected Band Structures (interpolated)
  - Density of States (DOS)
  - Screening parameters
  - Optimized Wannier orbitals

## Interfaces & Ecosystem
- **Quantum ESPRESSO**: Functions as a driver/wrapper around QE (`pw.x`, `ph.x`).
- **Wannier90**: Tightly coupled for orbital localization and interpolation.
- **ASE**: Compatible with Atomic Simulation Environment workflows.
- **Python API**: Provides a scriptable interface for advanced users.

## Performance Characteristics
- **Speed**: Significantly faster than GW (Many-Body Perturbation Theory). Roughly 2-10x cost of standard DFT depending on functional.
- **Scaling**: Linear-scaling characteristics arising from the localization of Wannier functions.
- **Parallelization**: Inherits MPI/OpenMP scalability from Quantum ESPRESSO engines.

## Limitations & Known Constraints
- **Metals**: Primarily designed for systems with a band gap (insulators, semiconductors, molecules); application to metals is non-trivial.
- **Dependence on Wannierization**: Quality of results depends on obtaining good Maximally Localized Wannier Functions (MLWFs).
- **Memory**: Can be memory-intensive for large supercells when calculating screening.

## Comparison with Other Codes
- **vs. GW (Yambo, BerkeleyGW)**: Koopmans is faster and avoids convergence issues with empty states/frequency integration, though GW is more general for metals.
- **vs. Hybrid Functionals (HSE06)**: Koopmans theory is non-empirical and often yields more accurate gaps than standard hybrids which rely on fixed mixing parameters.
- **vs. DFT+U**: DFT+U corrects local errors but doesn't guarantee correct spectral properties; Koopmans enforces the correct physical condition for charged excitations.

## Application Areas
- **Photovoltaics**: Accurate band alignment and gap prediction for solar cell materials.
- **Catalysis**: Correct placement of orbital levels in surface reactions.
- **Molecular Electronics**: Transport gaps in organic semiconductors.
- **Spectroscopy**: Interpretation of photoemission experiments.

## Community and Support
- **Development**: Led by the THEOS group (EPFL) and MARVEL NCCR.
- **Documentation**: Comprehensive ReadTheDocs with tutorials.
- **Forum**: Support via Quantum ESPRESSO forums and GitHub issues.

## Verification & Sources
- **Official Website**: [https://koopmans-functionals.org/](https://koopmans-functionals.org/)
- **Repository**: [https://github.com/koopmans-functionals/koopmans](https://github.com/koopmans-functionals/koopmans)
- **Primary Publication**: I. Dabo et al., Phys. Rev. B 82, 115121 (2010); N. Colonna based papers.
- **Verification status**: ✅ VERIFIED
  - Active GitHub repository.
  - Validated against experimental band gaps for large semiconductor datasets.
  - Part of the Quantum ESPRESSO ecosystem.
