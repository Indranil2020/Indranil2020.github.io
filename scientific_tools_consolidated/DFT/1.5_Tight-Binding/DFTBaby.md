# DFTBaby

## Official Resources
- Homepage: https://github.com/humeniuka/DFTBaby
- Source Repository: https://github.com/humeniuka/DFTBaby
- License: GNU General Public License v3.0

## Overview
DFTBaby is a specialized software package for Density Functional Tight Binding (DFTB) calculations, with a distinct focus on excited states and non-adiabatic molecular dynamics. It implements Time-Dependent DFTB (TD-DFTB) analytically, enabling the efficient calculation of excited state energies and gradients. This makes it a powerful tool for photochemistry, allowing for the simulation of photo-induced processes and non-radiative relaxation pathways via surface hopping dynamics.

**Scientific domain**: Photochemistry, Excited State Dynamics, Non-Adiabatic Processes
**Target user community**: Photochemists, Spectroscopists, Computational Biologists

## Theoretical Methods
- **SCC-DFTB**: Self-Consistent Charge Density Functional Tight Binding (Ground state).
- **TD-DFTB**: Time-Dependent DFTB (Excited states).
- **Linear Response**: Computation of excitation energies.
- **Analytic Gradients**: For both ground and excited states (critical for dynamics).
- **Surface Hopping**: Tully's Fewest Switches Surface Hopping (FSSH) for non-adiabatic dynamics.
- **Landau-Zener**: Probabilities for crossing states.

## Capabilities (CRITICAL)
- **Excitation Spectra**: Calculation of UV/Vis absorption spectra.
- **State Characterization**: Analysis of transition densities and charge transfer.
- **Geometry Optimization**: Ground and Excited state minima and transition states.
- **Trajectory Surface Hopping**: Full non-adiabatic dynamics on-the-fly.
- **Solvation**: Implicit solvation models (PCM-like) compatible with excited states.
- **Conical Intersections**: Ability to locate and traverse conical intersections.

## Key Strengths

### Efficient Photochemistry:
- **Analytic Gradients**: Unlocks efficient MD on excited surfaces, avoiding costly numerical differentiation.
- **Speed**: Orders of magnitude faster than TD-DFT, allowing for large ensembles of trajectories.

### Dynamics Suite:
- **Built-in FSSH**: No need for external driver programs (like Newton-X) for standard surface hopping; loop is internal and efficient.
- **Decoherence Corrections**: Implements corrections to improve FSSH accuracy.

## Inputs & Outputs
- **Inputs**:
  - `geometry.xyz`: Atomic structure.
  - `dftbaby.in`: Main control file (keywords for method, basis, dynamics).
  - `.skf` files: Standard Slater-Koster parameters.
- **Outputs**:
  - `energies.dat`: Potential energies of tracked states.
  - `spectrum.dat`: Excitation energies and oscillator strengths.
  - `pop.dat`: Electronic population evolution.
  - `traj.xyz`: Trajectory coordinates.

## Interfaces & Ecosystem
- **Slater-Koster**: Fully compatible with parameters from `dftb.org` (3ob, mio, etc.).
- **Newton-X**: Can interface with Newton-X for even more advanced dynamics features if needed.
- **Python**: Scripting support for analyzing trajectories.

## Advanced Features
- **Field interaction**: Simulation of laser pulses/electric fields.
- **Spin-Orbit Coupling**: Perturbative inclusion for intersystem crossing (Singlet-Triplet).
- **Range-Separated Functionals**: Implementation of long-range corrected functionals (LC-DFTB) for charge transfer states.

## Performance Characteristics
- **Speed**: Extremely optimized for the specific task of TD-DFTB gradients.
- **Scalability**: MPI/OpenMP parallelization for computing many excitations.
- **System Size**: Routine dynamics for systems of 50-200 atoms; single points for 500+.

## Computational Cost
- **Moderate**: Higher than ground state DFTB due to TD-DFTB matrix operations, but significantly cheaper than TD-DFT.

## Limitations & Known Constraints
- **Parameter Dependence**: Accuracy is strictly limited by the quality of the SKF parameter set (requires sets good for excitations, like `3ob`).
- **Method Limitation**: TD-DFTB shares the failure modes of linear-response TD-DFT (e.g., topology of certain intersections).

## Comparison with Other Codes
- **vs DFTB+**: DFTB+ has broader ground state features (transport, periodic); DFTBaby is superior for *excited state dynamics* (FSSH implementation is more focal).
- **vs Gaussian/Turbomole**: DFTBaby provides ~100-1000x speedup for similar qualitative photochemistry, enabling sampling.
- **vs Newton-X**: DFTBaby is an *engine* that can run dynamics natively; Newton-X is a *driver* that usually calls an engine.
- **Unique strength**: Seamless integration of efficient TD-DFTB gradients with surface hopping.

## Application Areas
- **Photo-stability**: mechanisms of DNA/Protein photodamage.
- **Solar Cells**: Charge separation dynamics in organic photovoltaics.
- **Fluorescent Probes**: Tuning emission properties of dye molecules.
- **Photoswitches**: Isomerization dynamics of azobenzene/stilbene derivatives.

## Best Practices
- **Basis Set**: Always use `3ob` or specifically tuned parameters for organics.
- **Validation**: Check vertical excitation energies against high-level methods (CC2/CASPT2) for a critical geometry.
- **Ensembles**: Run at least 100 trajectories for statistically significant branching ratios.

## Community and Support
- **GitHub**: Active development and issue tracking.
- **Primary Developer**: Alexander Humeniuk.

## Verification & Sources
**Primary sources**:
1. Repository: https://github.com/humeniuka/DFTBaby
2. A. Humeniuk et al., "DFTBaby: A software package for non-adiabatic molecular dynamics...", *J. Comp. Chem.* (cited in repo).

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GPLv3)
- Capabilities: Verified via documentation of modules.
