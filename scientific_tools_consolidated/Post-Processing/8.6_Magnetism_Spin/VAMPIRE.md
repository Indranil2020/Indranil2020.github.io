# VAMPIRE (Viscoelastic Atomistic Magnetic simulaton Program for Interactive Research and Engineering)

## Official Resources
- Homepage: https://vampire.york.ac.uk/
- Documentation: https://vampire.york.ac.uk/documentation/
- Source Repository: https://github.com/richard-evans/vampire
- License: GNU General Public License v3.0

## Overview
VAMPIRE is an atomistic spin dynamics (ASD) code designed for the simulation of magnetic materials. It allows for the simulation of magnetic properties at finite temperatures, including magnetization dynamics, hysteresis loops, and Curie temperatures. VAMPIRE uses the Landau-Lifshitz-Gilbert (LLG) equation and Monte Carlo methods to model magnetic systems ranging from single atoms to complex heterostructures.

**Scientific domain**: Atomistic spin dynamics, micromagnetics, magnetic materials  
**Target user community**: Magnetism researchers, recording media engineers, spintronics

## Theoretical Methods
- Landau-Lifshitz-Gilbert (LLG) equation
- Stochastic LLG (for finite temperature dynamics)
- Monte Carlo simulation (Metropolis)
- Constrained Monte Carlo
- Heisenberg Hamiltonian (Exchange, Anisotropy, Zeeman, Dipole-Dipole)
- Two-temperature model (for ultrafast demagnetization)

## Capabilities (CRITICAL)
- Calculation of Curie temperatures and temperature-dependent magnetization
- Simulation of hysteresis loops and coercive fields
- Ultrafast demagnetization dynamics (laser heating)
- Exchange bias and interface effects
- Core-shell nanoparticles and heterostructures
- High-performance parallelization (MPI)
- Analysis of spin configurations

**Sources**: VAMPIRE website, J. Phys.: Condens. Matter 26, 103202 (2014)

## Inputs & Outputs
- **Input formats**: `vampire.input` (control parameters), `vampire.mat` (material properties), `vampire.UCF` (unit cell)
- **Output data types**: `output` (magnetization vs time/field), `atoms` (spin snapshots), `susceptibility`

## Interfaces & Ecosystem
- **Visualization**: Output compatible with POV-Ray, JMOL, RasMol
- **TB2J**: Can use exchange parameters calculated by TB2J
- **Scripts**: Python/Bash scripts for managing simulations

## Workflow and Usage
1. Define crystal structure (UCF file).
2. Define material properties (exchange J, anisotropy K).
3. Configure simulation type in `vampire.input` (e.g., Curie temperature, Hysteresis).
4. Run `vampire-serial` or `vampire-parallel`.
5. Analyze output data files.

## Performance Characteristics
- Highly optimized C++ code
- Parallelized with MPI for large systems (millions of spins)
- Efficient neighbor lists

## Application Areas
- Magnetic recording media (HAMR)
- Permanent magnets
- Magnetic nanoparticles (hyperthermia)
- Spintronic devices
- Antiferromagnets

## Community and Support
- Open-source (GPL v3)
- Developed at University of York (Richard Evans group)
- Active development and user support

## Verification & Sources
**Primary sources**:
1. Homepage: https://vampire.york.ac.uk/
2. GitHub: https://github.com/richard-evans/vampire
3. Publication: R. F. L. Evans et al., J. Phys.: Condens. Matter 26, 103202 (2014)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN (GitHub)
- Development: ACTIVE (Univ. York)
- Applications: Spin dynamics, LLG, magnetic materials, recording media
