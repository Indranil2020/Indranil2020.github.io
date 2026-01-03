# XtalOpt

## Official Resources
- Homepage: https://xtalopt.github.io/
- Documentation: https://xtalopt.github.io/documentation.html
- Source Repository: https://github.com/xtalopt/XtalOpt
- License: GNU General Public License v3.0

## Overview
XtalOpt is an open-source evolutionary algorithm for crystal structure prediction. It is implemented as an extension to the Avogadro molecular editor, providing a graphical user interface for setting up and monitoring structure prediction runs. XtalOpt searches for stable crystal structures by optimizing randomly generated structures and evolving them through genetic operations like crossover and mutation.

**Scientific domain**: Crystal structure prediction, evolutionary algorithms  
**Target user community**: Materials scientists, chemists, crystallographers

## Theoretical Methods
- Evolutionary Algorithm (EA)
- Genetic operations: Crossover, Mutation, Strain
- Local optimization via external codes
- Random initialization with space group constraints
- Duplicate removal via fingerprinting

## Capabilities (CRITICAL)
- Crystal structure prediction from composition
- User-friendly GUI within Avogadro
- Real-time visualization of the search progress
- Interactive plot of energy vs structure
- Support for various unit cell constraints
- Duplicate structure detection (XtalComp)
- Integration with multiple optimizers (VASP, GULP, PWscf, SIESTA, CASTEP, ATK, etc.)

**Sources**: XtalOpt website, Comp. Phys. Comm. 182, 372 (2011)

## Inputs & Outputs
- **Input formats**: GUI-based setup, optimizer templates
- **Output data types**: Optimized structures (.cif, .xyz), search history, energy plots

## Interfaces & Ecosystem
- **Avogadro**: Host application for XtalOpt
- **Optimizers**: VASP, GULP, Quantum ESPRESSO, etc.
- **XtalComp**: Structure comparison library

## Workflow and Usage
1. Open Avogadro and start XtalOpt extension.
2. Define composition and simulation parameters (population size, number of generations).
3. Configure optimization steps (code, potentials).
4. Start search.
5. Monitor progress in GUI: visualize structures, check convergence.

## Performance Characteristics
- Efficiency depends on the underlying optimizer.
- GUI overhead is minimal.
- Parallelization managed by submitting jobs to clusters (PBS/Torque/SLURM support).

## Application Areas
- New material discovery
- High-pressure phases
- Hydrides and oxides
- Superconductors

## Community and Support
- Open-source (GPL v3)
- GitHub repository
- Developed by Zurek Group (SUNY Buffalo)

## Verification & Sources
**Primary sources**:
1. Homepage: https://xtalopt.github.io/
2. GitHub: https://github.com/xtalopt/XtalOpt
3. Publication: D.C. Lonie and E. Zurek, Comp. Phys. Comm. 182, 372 (2011)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: AVAILABLE
- Source: OPEN (GitHub)
- Development: ACTIVE (Zurek Group)
- Applications: Evolutionary structure prediction, GUI-based, open source
