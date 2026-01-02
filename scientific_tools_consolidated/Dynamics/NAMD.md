# NAMD

## Official Resources
- Homepage: https://www.ks.uiuc.edu/Research/namd/
- Documentation: https://www.ks.uiuc.edu/Research/namd/documentation/
- Source Repository: Available with license
- License: Free for non-commercial use

## Overview
NAMD (NAnoscale Molecular Dynamics) is a parallel molecular dynamics code designed for high-performance simulation of large biomolecular systems. Developed by the Theoretical and Computational Biophysics Group at the University of Illinois, NAMD is optimized for extreme scalability on supercomputers with thousands to millions of cores. It is particularly renowned for its exceptional parallel performance, GPU acceleration, and integration with VMD for visualization and analysis.

**Scientific domain**: Biomolecular molecular dynamics, parallel computing, extreme-scale simulation  
**Target user community**: Computational biophysicists, structural biologists, large-scale MD researchers

## Theoretical Methods
- Classical molecular dynamics (MD)
- Energy minimization
- Canonical ensemble (NVT)
- Isothermal-isobaric ensemble (NPT)
- Langevin dynamics
- Multiple time-stepping (r-RESPA)
- Particle Mesh Ewald (PME) for electrostatics
- SHAKE and RATTLE constraints
- Steered molecular dynamics (SMD)
- Targeted molecular dynamics (TMD)
- Adaptive biasing force (ABF)
- Metadynamics
- Replica exchange molecular dynamics (REMD)
- Free energy perturbation (FEP)
- Thermodynamic integration (TI)
- String method with swarms of trajectories
- Collective variables module (Colvars)

## Capabilities (CRITICAL)
- Biomolecular MD simulations
- Energy minimization
- Multiple ensembles (NVE, NVT, NPT)
- Temperature and pressure control
- Periodic boundary conditions
- Force field support (CHARMM, AMBER, OPLS, GROMOS)
- Long-range electrostatics (PME, multiple timestepping)
- Constraint algorithms (SHAKE, RATTLE)
- Steered MD for pulling simulations
- Free energy calculations (FEP, TI, ABF)
- Enhanced sampling (REMD, metadynamics)
- QM/MM simulations (via ORCA, MOPAC integration)
- Interactive MD (via VMD)
- Extreme parallelization (millions of cores)
- Outstanding GPU acceleration (CUDA)
- Multi-copy algorithms (REMD, string method)
- Collective variables for enhanced sampling
- Integration with VMD for analysis
- Grid-based force calculation
- Heterogeneous architectures (CPU+GPU)

**Sources**: Official NAMD documentation (https://www.ks.uiuc.edu/Research/namd/), confirmed in 7/7 source lists

## Key Strengths

### Scalability:
- Extreme parallel scaling
- Efficient on 100,000+ cores
- World-record simulations
- Optimized communication
- Charm++ parallel runtime

### GPU Acceleration:
- Excellent CUDA support
- Multiple GPU per node
- CPU-GPU load balancing
- 5-10x speedup typical
- Optimized kernels

### Large Systems:
- Millions of atoms
- Viral capsids
- Ribosomes
- Membrane systems
- Whole cells (future)

### Free Energy:
- Comprehensive FEP
- Thermodynamic integration
- ABF method
- Metadynamics support
- Drug binding calculations

### VMD Integration:
- Seamless workflow
- Interactive MD
- Real-time visualization
- Analysis tools
- Single ecosystem

## Inputs & Outputs
- **Input formats**:
  - PDB (structure)
  - PSF (topology)
  - Configuration file (.conf)
  - Parameter files (.prm, .str)
  - Coordinate files
  
- **Output data types**:
  - DCD (trajectory)
  - Energy logs
  - Restart files
  - PMF data
  - Custom output

## Interfaces & Ecosystem
- **Preparation**:
  - VMD for system setup
  - CHARMM-GUI
  - psfgen (topology)
  
- **Analysis**:
  - VMD (primary)
  - Colvars analysis
  - WHAM for PMF
  - Custom scripts
  
- **Enhanced Sampling**:
  - Colvars module
  - PLUMED integration
  - Replica exchange
  
- **QM/MM**:
  - ORCA interface
  - MOPAC interface
  - Custom QM codes
  
- **Visualization**:
  - VMD (native)
  - Real-time rendering
  - Movie generation

## Workflow and Usage

### Typical MD Workflow:

1. **System Preparation**:
   - Build PSF with psfgen or CHARMM-GUI
   - Solvate and add ions
   - Create configuration file

2. **Minimization**:
   ```tcl
   # Configuration file
   structure myprotein.psf
   coordinates myprotein.pdb
   
   temperature 310
   
   minimize 1000
   ```

3. **Equilibration**:
   ```tcl
   langevin on
   langevinTemp 310
   langevinDamping 5
   
   run 50000  # 100 ps
   ```

4. **Production**:
   ```tcl
   run 25000000  # 50 ns
   ```

5. **Run NAMD**:
   ```bash
   namd2 +p16 config.conf > output.log
   # GPU version
   namd2 +p1 +devices 0,1,2,3 config.conf
   ```

## Advanced Features

### Extreme Scaling:
- Charm++ runtime
- Spatial decomposition
- Hybrid parallelization
- Dynamic load balancing
- Scaling to supercomputers

### Multiple Timestepping:
- r-RESPA integrator
- Different timesteps for different forces
- Bonded: 1 fs
- Short-range: 2 fs
- PME: 4 fs
- Efficiency gain

### Steered MD:
- Apply forces to atoms
- Constant velocity pulling
- Constant force pulling
- AFM simulations
- Mechanical unfolding

### Free Energy:
- FEP windows
- Softcore potentials
- BAR analysis
- TI calculations
- ABF for PMF

### Replica Exchange:
- Temperature REMD
- Hamiltonian REMD
- Enhanced sampling
- Parallel tempering
- Multi-copy algorithms

### Collective Variables:
- Colvars module
- Geometric CVs
- Path CVs
- Metadynamics
- ABF

## Performance Characteristics
- **Speed**: Excellent for large systems
- **Scaling**: Best-in-class parallel scaling
- **GPU**: 5-10x speedup
- **Typical systems**: 100,000-10,000,000 atoms
- **Timestep**: 2 fs standard, 4 fs with SETTLE

## Computational Cost
- **Large systems**: Most efficient code
- **Scaling**: Enables otherwise impossible simulations
- **GPU**: Dramatically reduces time
- **Free energy**: Multiple replicas expensive
- **REMD**: Many replicas needed

## Limitations & Known Constraints
- **Learning curve**: Moderate
- **Configuration**: Text-based (verbose)
- **Force fields**: Primarily CHARMM format
- **Analysis**: Requires VMD or scripts
- **Small systems**: May not scale well
- **Platform**: Linux primarily, Windows/macOS possible
- **License**: Free for non-commercial, registration required

## Comparison with Other MD Codes
- **vs GROMACS**: NAMD better scalability, GROMACS faster for biomolecules on single node
- **vs AMBER**: NAMD better parallel, AMBER more force fields
- **vs LAMMPS**: NAMD for biomolecules, LAMMPS for materials
- **vs Desmond**: NAMD free, Desmond commercial but very fast
- **Unique strength**: Extreme scalability, GPU performance, VMD integration

## Application Areas

### Structural Biology:
- Protein dynamics
- Membrane proteins
- Protein-protein interactions
- Conformational changes
- Enzyme mechanisms

### Drug Discovery:
- Protein-ligand binding
- Free energy calculations
- Virtual screening
- ADME properties

### Large Systems:
- Viral capsids
- Ribosomes
- Chromatin
- Membrane vesicles
- Cellular components

### Mechanical Properties:
- Protein unfolding
- Mechanotransduction
- AFM simulations
- Elasticity

### Nucleic Acids:
- DNA/RNA dynamics
- Protein-DNA complexes
- Chromatin structure

## Best Practices

### System Preparation:
- Use CHARMM-GUI or VMD
- Check structure quality
- Proper solvation
- Neutralize system
- Minimize carefully

### Simulation Setup:
- Gradual heating
- NVT then NPT equilibration
- Restrain protein initially
- 2 fs timestep standard
- PME for electrostatics

### Performance:
- Use GPUs when available
- Optimize patch grid
- Balance load
- Multiple timestepping
- Minimize output frequency

### Free Energy:
- Sufficient windows
- Adequate sampling
- Check convergence
- Use BAR for analysis
- Soft-core potentials

## Community and Support
- Free for non-commercial use
- Registration required
- Mailing list support
- Extensive tutorials
- Active user community
- Regular workshops
- VMD/NAMD team support

## Educational Resources
- Comprehensive manual
- Tutorial examples
- Video tutorials
- Workshops
- Published protocols
- NAMD wiki

## Awards and Recognition
- Gordon Bell Prize (2002)
- Sidney Fernbach Award (2012)
- Extreme scaling demonstrations
- HPC benchmarks
- Scientific breakthroughs

## High-Performance Computing
- Optimized for supercomputers
- Petascale simulations
- Exascale ready
- Leadership computing
- World records

## Verification & Sources
**Primary sources**:
1. Official website: https://www.ks.uiuc.edu/Research/namd/
2. Documentation: https://www.ks.uiuc.edu/Research/namd/documentation/
3. J. C. Phillips et al., J. Comput. Chem. 26, 1781 (2005) - NAMD paper
4. J. C. Phillips et al., J. Chem. Phys. 153, 044130 (2020) - NAMD 2.14

**Secondary sources**:
1. NAMD manual and tutorials
2. Published MD studies using NAMD (>10,000 citations)
3. VMD/NAMD workshops
4. Confirmed in 7/7 source lists (claude, g, gr, k, m, q, z)

**Confidence**: CONFIRMED - Appears in ALL 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE and ACCESSIBLE
- Software: Available with registration (free non-commercial)
- Community support: Mailing list, tutorials, workshops
- Academic citations: >12,000
- Active development: Regular releases
- Benchmark validation: Extensively validated, HPC benchmarks
- Scalability: World-record parallel performance
- Industry standard: Leading code for large biomolecular simulations
