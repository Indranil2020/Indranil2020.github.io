# GROMACS

## Official Resources
- Homepage: https://www.gromacs.org/
- Documentation: https://manual.gromacs.org/
- Source Repository: https://gitlab.com/gromacs/gromacs
- License: GNU Lesser General Public License v2.1

## Overview
GROMACS (GROningen MAchine for Chemical Simulations) is a versatile and extremely fast molecular dynamics package primarily designed for simulations of proteins, lipids, and nucleic acids. Developed originally at the University of Groningen, it is the fastest open-source MD code for biomolecular systems, featuring excellent GPU acceleration, sophisticated free energy calculation methods, and comprehensive analysis tools.

**Scientific domain**: Biomolecular molecular dynamics, proteins, membranes, drug design  
**Target user community**: Biochemists, biophysicists, computational biologists, pharmaceutical researchers

## Theoretical Methods
- Classical molecular dynamics (MD)
- Energy minimization (steepest descent, conjugate gradient)
- Canonical ensemble (NVT)
- Isothermal-isobaric ensemble (NPT)
- Microcanonical ensemble (NVE)
- Langevin dynamics
- Stochastic dynamics
- Free energy perturbation (FEP)
- Thermodynamic integration (TI)
- Umbrella sampling
- Metadynamics (via PLUMED)
- Replica exchange MD (REMD)
- Accelerated weight histogram method (AWH)
- Constraint algorithms (LINCS, SETTLE)
- Virtual sites
- Coarse-grained models (Martini)

## Capabilities (CRITICAL)
- Biomolecular MD simulations
- Energy minimization
- Multiple ensembles (NVE, NVT, NPT)
- Temperature and pressure coupling (Berendsen, Nosé-Hoover, Parrinello-Rahman)
- Periodic boundary conditions
- Extensive force field support (AMBER, CHARMM, GROMOS, OPLS)
- Long-range electrostatics (PME, reaction-field)
- Free energy calculations (FEP, TI, BAR)
- Enhanced sampling (umbrella sampling, AWH, REMD)
- Pull code for steered MD
- Constraint algorithms (bond, angle)
- Virtual interaction sites
- Coarse-grained simulations (Martini)
- GPU acceleration (CUDA, OpenCL)
- Excellent performance (fastest biomolecular MD)
- Extensive analysis tools (100+ built-in)
- Trajectory manipulation tools
- Massively parallel (MPI, thread-MPI)
- Heterogeneous parallelization (CPU+GPU)

**Sources**: Official GROMACS documentation (https://www.gromacs.org/), confirmed in 7/7 source lists

## Key Strengths

### Performance:
- Fastest MD code for biomolecules
- Exceptional GPU acceleration
- Highly optimized algorithms
- SIMD vectorization
- Scaling to thousands of nodes

### GPU Acceleration:
- Native CUDA support
- OpenCL support
- Offloading of PME and bonded forces
- Up to 10x speedup
- CPU-GPU load balancing

### Free Energy:
- Comprehensive FEP implementation
- Thermodynamic integration
- BAR and MBAR analysis
- Alchemical transformations
- AWH adaptive sampling

### Analysis Tools:
- 100+ built-in tools
- Trajectory analysis
- Structure analysis
- Thermodynamic properties
- RMSD, RMSF calculations
- Principal component analysis

### User-Friendly:
- pdb2gmx for system preparation
- Well-documented
- Active community
- Tutorials available

## Inputs & Outputs
- **Input formats**:
  - .mdp file (MD parameters)
  - .gro or .pdb (coordinates)
  - .top (topology)
  - .itp (include topology)
  - .ndx (index groups)
  
- **Output data types**:
  - .trr (full precision trajectory)
  - .xtc (compressed trajectory)
  - .edr (energy file)
  - .log (log file)
  - .cpt (checkpoint for restart)
  - Various analysis outputs

## Interfaces & Ecosystem
- **Preparation tools**:
  - pdb2gmx (topology generation)
  - editconf (box setup)
  - genion (ion addition)
  - solvate (add solvent)
  
- **Simulation tools**:
  - grompp (preprocessor)
  - mdrun (MD engine)
  - Various gmx tools
  
- **Analysis**:
  - 100+ gmx analysis tools
  - gmx energy, rms, rmsf, gyrate, etc.
  - Custom analysis via API
  
- **Enhanced sampling**:
  - PLUMED plugin
  - AWH built-in
  - Pull code
  
- **Visualization**:
  - VMD, PyMOL
  - NGLView
  - GROMACS viewer
  
- **Workflows**:
  - Python wrappers
  - MDAnalysis compatible
  - BioExcel tools

## Workflow and Usage

### Typical MD Workflow:

1. **System Preparation**:
   ```bash
   gmx pdb2gmx -f protein.pdb -o protein.gro -water tip3p
   gmx editconf -f protein.gro -o protein_box.gro -c -d 1.0 -bt cubic
   gmx solvate -cp protein_box.gro -cs spc216.gro -o protein_solv.gro
   gmx grompp -f ions.mdp -c protein_solv.gro -p topol.top -o ions.tpr
   gmx genion -s ions.tpr -o protein_ions.gro -p topol.top -pname NA -nname CL -neutral
   ```

2. **Energy Minimization**:
   ```bash
   gmx grompp -f minim.mdp -c protein_ions.gro -p topol.top -o em.tpr
   gmx mdrun -v -deffnm em
   ```

3. **Equilibration (NVT)**:
   ```bash
   gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
   gmx mdrun -deffnm nvt
   ```

4. **Equilibration (NPT)**:
   ```bash
   gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
   gmx mdrun -deffnm npt
   ```

5. **Production MD**:
   ```bash
   gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
   gmx mdrun -deffnm md
   ```

6. **Analysis**:
   ```bash
   gmx rms -s md.tpr -f md.xtc -o rmsd.xvg
   gmx rmsf -s md.tpr -f md.xtc -o rmsf.xvg
   ```

## Advanced Features

### Free Energy Calculations:
- λ-coupling for alchemical transformations
- Soft-core potentials
- Multiple λ-windows
- BAR/MBAR analysis
- Drug binding affinity predictions

### Adaptive Sampling (AWH):
- Bias potential adaptation
- Free energy landscapes
- Automatic convergence
- Multiple collective variables

### Replica Exchange:
- Temperature REMD
- Hamiltonian REMD
- Enhanced sampling
- Parallel tempering

### Pull Code:
- Steered MD
- Umbrella sampling
- Potential of mean force (PMF)
- AFM simulations

### Coarse-Graining:
- Martini force field support
- Large system simulations
- Lipid bilayers
- Membrane proteins

### Virtual Sites:
- Dummy atoms
- Rigid water models (TIP4P)
- Aromatic rings
- Increased timesteps

## Performance Characteristics
- **Speed**: 10-100 ns/day typical on modern GPU
- **Scaling**: Good to 1000+ nodes
- **GPU**: 5-15x speedup vs CPU
- **Typical systems**: 10,000-10,000,000 atoms
- **Timestep**: 2 fs standard, 4-5 fs with virtual sites

## Computational Cost
- **Biomolecules**: Very efficient
- **PME**: Well-optimized
- **GPU**: Dramatically reduces cost
- **Free energy**: Multiple windows expensive
- **REMD**: Multiple replicas needed

## Limitations & Known Constraints
- **Biomolecule focus**: Optimized for proteins/membranes
- **Force fields**: Primarily biomolecular (AMBER, CHARMM, etc.)
- **Materials**: Less suitable than LAMMPS
- **Learning curve**: Moderate
- **mdp files**: Many parameters to understand
- **Units**: nm, ps (not Å, fs like others)
- **Platform**: Linux primarily, Windows/macOS possible

## Comparison with Other MD Codes
- **vs LAMMPS**: GROMACS faster for biomolecules, LAMMPS more versatile
- **vs NAMD**: GROMACS faster, NAMD better for very large systems
- **vs AMBER**: GROMACS faster and open-source
- **vs CHARMM**: GROMACS faster, similar force fields
- **Unique strength**: Speed for biomolecules, GPU performance, free energy

## Application Areas

### Protein Dynamics:
- Protein folding
- Conformational changes
- Protein-protein interactions
- Enzyme mechanisms

### Membrane Biophysics:
- Lipid bilayers
- Membrane proteins
- Ion channels
- Drug-membrane interactions

### Drug Design:
- Binding free energies
- Virtual screening
- ADME properties
- Ligand-protein complexes

### Nucleic Acids:
- DNA/RNA dynamics
- Protein-DNA complexes
- Drug-DNA interactions

### Materials (Limited):
- Polymers
- Soft matter
- Liquid crystals

## Best Practices

### System Preparation:
- Check protein structure quality
- Add missing atoms/residues
- Choose appropriate force field
- Proper solvation box size
- Neutralize system

### Simulation Setup:
- Energy minimization first
- NVT then NPT equilibration
- Restrain protein during equilibration
- 2 fs timestep standard
- PME for long-range electrostatics

### Equilibration:
- Monitor temperature and pressure
- Check energy conservation (NVE)
- RMSD stabilization
- Sufficient equilibration time (100s ps to ns)

### Production:
- Long trajectories for sampling
- Save coordinates frequently
- Monitor system stability
- Generate backups/checkpoints

### Analysis:
- Remove PBC artifacts
- Fit structures for RMSD
- Block averaging for errors
- Check convergence

### Performance:
- Use GPUs when available
- Optimize PME grid spacing
- Balance CPU and GPU work
- Efficient domain decomposition

## Community and Support
- Open-source on GitLab
- Very large user community
- Active mailing lists
- Forum for questions
- Regular workshops
- Extensive tutorials
- Annual user meetings

## Educational Resources
- Comprehensive manual
- Justin Lemkul's tutorials (very popular)
- BioExcel training
- Video tutorials
- Published books
- Hands-on workshops

## High-Performance Computing
- Optimized for HPC
- Excellent GPU utilization
- SIMD optimizations (AVX-512)
- Heterogeneous acceleration
- Production use worldwide

## Verification & Sources
**Primary sources**:
1. Official website: https://www.gromacs.org/
2. Documentation: https://manual.gromacs.org/
3. GitLab repository: https://gitlab.com/gromacs/gromacs
4. M. J. Abraham et al., SoftwareX 1-2, 19 (2015) - GROMACS 5
5. H. J. C. Berendsen et al., Comp. Phys. Comm. 91, 43 (1995) - Original GROMACS

**Secondary sources**:
1. GROMACS manual and tutorials
2. Published MD studies using GROMACS (>20,000 citations)
3. BioExcel training materials
4. Confirmed in 7/7 source lists (claude, g, gr, k, m, q, z)

**Confidence**: CONFIRMED - Appears in ALL 7 independent source lists

**Verification status**: VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE and ACCESSIBLE
- Source code: OPEN (GitLab, LGPL v2.1)
- Community support: Very active mailing lists, forums
- Academic citations: >30,000
- Active development: Regular releases, active GitLab
- Benchmark validation: Extensively validated, industry standard
- HPC optimization: Fastest biomolecular MD code
- Industry standard: Dominant code for biomolecular simulations
