# McPhase

## Official Resources
- Homepage: https://www2.cpfs.mpg.de/~rotter/homepage_mcphase/
- Documentation: https://www2.cpfs.mpg.de/~rotter/homepage_mcphase/manual/manual.html
- Source Repository: Distributed via website (GPL)
- License: GNU General Public License v3.0

## Overview
McPhase is a software package for calculating magnetic phase diagrams and thermodynamic properties of magnetic materials. It uses a mean-field approximation combined with Monte Carlo simulations to treat localized magnetic moments in crystal fields. It is particularly specialized for rare-earth magnetism, handling crystal electric field (CEF) effects, exchange interactions, and magneto-elastic coupling.

**Scientific domain**: Magnetic phase diagrams, rare-earth magnetism, crystal fields  
**Target user community**: Solid-state physicists, neutron scattering researchers

## Theoretical Methods
- Mean Field Theory (for magnetic ordering)
- Crystal Electric Field (CEF) theory (Stevens operators)
- Monte Carlo simulation (for phase diagrams)
- Random Phase Approximation (RPA) for excitations
- Dipole approximation for magnetic moments
- Multipole interactions

## Capabilities (CRITICAL)
- Calculation of magnetic phase diagrams (H-T, P-T)
- Simulation of specific heat, magnetization, susceptibility
- Calculation of magnetic excitations (spin waves, crystal field levels)
- Neutron scattering cross-sections (elastic and inelastic)
- Handling of complex magnetic structures (commensurate/incommensurate)
- Magnetostriction and magneto-elastic effects

**Sources**: McPhase website, J. Magn. Magn. Mater. 272, E481 (2004)

## Inputs & Outputs
- **Input formats**: `.mcphase` (control), `.sipf` (interaction parameters), `.j1` (exchange constants)
- **Output data types**: `.out` (thermodynamics), `.calc` (calculated properties), `.jpg/png` (phase diagrams)

## Interfaces & Ecosystem
- **Simplot**: Built-in plotting tool
- **Display**: Visualization of magnetic structures
- **Mcdisp**: Module for dispersion relations
- **Fit**: Fitting module for experimental data

## Workflow and Usage
1. Define crystal structure and magnetic ions.
2. Parameterize Crystal Field (point charge or DFT fitted).
3. Define exchange interactions.
4. Run `mcphas` to calculate free energy and stabilize magnetic structure.
5. Run `mcdisp` for excitations.

## Performance Characteristics
- Mean field calculations are fast
- Monte Carlo steps can be time-consuming for large phase diagrams
- Efficient for localized moment systems

## Application Areas
- Rare-earth intermetallics
- Multiferroics
- Frustrated magnets
- Interpretation of neutron scattering data
- Magnetocaloric effect

## Community and Support
- Developed by Martin Rotter (Dresden/McPhase group)
- Academic user base
- Comprehensive manual and examples

## Verification & Sources
**Primary sources**:
1. Homepage: https://www2.cpfs.mpg.de/~rotter/homepage_mcphase/
2. Publication: M. Rotter et al., J. Magn. Magn. Mater. 272, E481 (2004)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: AVAILABLE
- Source: OPEN (GPL)
- Development: ACTIVE (Rotter)
- Applications: Magnetic phase diagrams, CEF, rare-earths
