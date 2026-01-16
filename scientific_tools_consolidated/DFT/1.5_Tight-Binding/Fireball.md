# Fireball

## Official Resources
- Homepage: https://github.com/FIREBALL2020
- Documentation: https://thunder-dft.github.io/
- Source Repository: https://github.com/FIREBALL2020/thunder-master
- License: GPLv3

## Overview
Fireball is an efficient *ab initio* tight-binding Density Functional Theory (DFT) code designed for molecular dynamics simulations of large systems. It utilizes a local-orbital formulation of DFT, enabling the simulation of supercells containing thousands of atoms with linear-scaling O(N) computational cost.

**Scientific domain**: Materials science, surface science, biomolecules, nanostructures
**Target user community**: Researchers performing large-scale MD simulations and studying extended systems

## Theoretical Methods
- Density Functional Theory (DFT)
- Tight-Binding (TB) formalism
- Local-orbital basis sets (Numerical Atomic Orbitals)
- Sankey-Niklewski (SN) approach
- Pseudopotentials (Norm-conserving)
- Molecular Dynamics (MD)
- Lewis-structure initialization

## Capabilities (CRITICAL)
- **Large-scale simulations**: Efficiently handles thousands of atoms (up to 10,000+).
- **Linear scaling**: O(N) cost for total energy and forces, avoiding cubic scaling of standard DFT.
- **Basis sets**: Optimized numerical atomic orbitals (NAOs), "Fireballs", which are strictly localized.
- **Dynamical reconstruction**: Suitable for surface catalytic processes and reconstruction.
- **Charge transport**: Investigations in amorphous systems (e.g., DNA) and molecular junctions.
- **QM/MM**: Interface with AMBER for hybrid quantum/classical simulations of biomolecules.
- **Van der Waals corrections**: Implementation available in newer versions (Fireball2020).
- **Electronic structure**: Band structures, Density of States (DOS).
- **Transport**: Conductance and STM image simulations.

**Sources**: Official GitHub (https://github.com/FIREBALL2020), cited in 4 sources.

## Key Strengths

### Efficiency via Local Orbitals:
- **Fireball Orbitals**: Numerical atomic orbitals with strict cutoff radii.
- **Pre-computed Integrals**: Three-center integrals are pre-calculated and stored, speeding up runtime.
- **Sparse Hamiltonian**: Local nature leads to sparse matrices, enabling O(N) scaling.

### Large-Scale MD:
- **Mesoscale Systems**: Designed for systems too large for plane-wave DFT but requiring quantum accuracy.
- **Long Timescales**: Efficiency allows for longer molecular dynamics trajectories.
- **Dynamical Evolution**: Ideal for studying kinetic processes and surface restructuring.

### Hybrid Applications (QM/MM):
- **Biomolecular Focus**: Strong integration with AMBER for enzymatic reactions and DNA studies.
- **Active Site Treatment**: Quantum treatment of active sites with classical environment.

## Inputs & Outputs
- **Input formats**:
  - `fireball.in`: Main control file (SCF settings, time steps, temperature).
  - `structure.inp`: Coordinate file.
  - Basis set files (pre-generated).
  - Pseudopotential files (Fdata).
  
- **Output data types**:
  - `param.dat`: Output parameters.
  - `dynamics.dat`: MD trajectory.
  - Energies and forces logs.
  - Electronic structure data (eigenvalues, DOS).
  - STM images (if requested).

## Interfaces & Ecosystem
- **QM/MM Integration**:
  - **AMBER**: Direct interface for hybrid calculations.
  
- **Tools**:
  - **Lightning**: Fast visualization tool/pre-processor.
  - **ASE (Atomic Simulation Environment)**: Scripting and workflow control.
  - **FireballTG**: Fork with enhanced transport capabilities (Transiesta-like).

## Workflow and Usage

### Typical MD Workflow:
1. **Preparation**: Generate initial structure, select basis set and pseudopotentials.
2. **Setup**: Configure `fireball.in` for `iensemble` (NVE/NVT) and `dt` (timestep).
3. **Execution**: Run Fireball (parallel execution supported).
4. **Analysis**: Extract `dynamics.dat` for visualization and structural analysis.

### Transport Workflow:
1. **Lead Definition**: Define semi-infinite leads and scattering region.
2. **Calculation**: Compute Green's functions.
3. **Output**: Transmission spectra and I-V curves.

## Advanced Features

### Surface Science:
- **Dynamical Reconstruction**: Observing surface atom rearrangement in real-time.
- **STM Simulation**: Generating theoretical Scanning Tunneling Microscopy images.
- **NEB**: Nudged Elastic Band for reaction barriers (in some versions).

### Transport (FireballTG):
- **NEGF Formalism**: Non-Equilibrium Green's Function for molecular electronics.
- **Conductance**: Landauer-Buttiker formalism.

## Performance Characteristics
- **Speed**: Significantly faster than plane-wave DFT (e.g., VASP, QE) for large systems.
- **Scaling**: Strictly linear O(N) for matrix construction; diagonalization can be O(N^3) or O(N) depending on solver.
- **System Size**: Routine handling of 1,000-5,000 atoms.

## Computational Cost
- **Single-Point**: Seconds to minutes for hundreds of atoms.
- **MD Step**: Rapid enough for ps-scale dynamics in reasonable wall time.
- **Comparison**: Slower than DFTB+ (which uses parameter tables), but more rigorous (explicit integrals).

## Limitations & Known Constraints
- **Accuracy vs. Efficiency**: "Fireball" approximations (basis cutoff) can lead to reduced accuracy compared to converged plane-wave results.
- **Self-Consistency**: Some older versions or modes (non-SCC) lack full charge self-consistency, though modern versions address this.
- **Basis Set Optimization**: Requires careful selection of orbital radii; poor choices lead to errors.
- **TDDFT**: Limited excited state capabilities compared to specialized codes.
- **Documentation**: Can be fragmented between different versions (Fireball vs Fireball2020 vs FireballTG).

## Comparison with Other Codes
- **vs DFTB+**: Fireball calculates integrals *ab initio* (using efficient numerics) rather than using fitted tables. It is generally more transferable but computationally heavier than parameterized DFTB.
- **vs SIESTA**: Both use numerical atomic orbitals. SIESTA is more general-purpose; Fireball is highly specialized for MD speed in large systems.
- **vs VASP/QE**: Fireball is much faster for >500 atoms, but less accurate for high-precision electronic structure (e.g., band gaps).

## Application Areas

### Biomolecular Simulations:
- **Enzymatic Catalysis**: QM/MM study of reaction mechanisms.
- **DNA Damage**: Photolesions and charge transfer in DNA.
- **Protein Solvation**: Water-protein semi-quantum interactions.

### Nanomaterials:
- **Surface Reconstruction**: Semiconductor surface annealing.
- **Nanowires**: Conductance and structural stability.
- **Defects**: Diffusion of defects in bulk materials.

### Molecular Electronics:
- **Single-Molecule Junctions**: I-V characteristics of organic molecules between leads.

## Best Practices

### Basis Set Selection:
- **Validation**: Test numerical basis sets against plane-wave results for small systems first.
- **Cutoffs**: Ensure orbital cutoff radii are sufficient to capture bonding but short enough for efficiency.

### Convergence:
- **SCF**: Use mixing (Pulay/Broyden) cautiously; difficult systems may require increased electronic temperature (smearing).
- **MD Stability**: Use smaller timesteps (0.5-1.0 fs) for variable-charge/flexible-basis MD to conserve energy.

### QM/MM Setup:
- **Partitioning**: carefully define the QM region to minimize boundary errors.
- **Link Atoms**: Use standard link atom schemes (hydrogen caps) at the QM/MM boundary.

## Community and Support
- **Project Structure**: Open-source, hosted on GitHub.
- **Primary Group**: Fireball2020 / Thunder-DFT team.
- **Support**: GitHub issues, academic collaboration networks.

## Verification & Sources
**Primary sources**:
1. Official GitHub: https://github.com/FIREBALL2020
2. Documentation: https://thunder-dft.github.io/
3. Comparative Reviews: Included in lists of O(N) DFT codes.

**Confidence**: VERIFIED
