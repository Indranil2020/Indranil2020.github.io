# EXESS

## Official Resources
- Homepage: https://barcagrp.com/exess/
- Documentation: https://barcagrp.com/exess/
- Source Repository: Closed (Waitlist/Academic access)
- License: Academic/Commercial

## Overview
EXESS (Extreme-scale Electronic Structure System) is a GPU-native quantum chemistry code designed for extreme-scale ab initio molecular dynamics (AIMD) capabilities. It won the 2024 ACM Gordon Bell Prize for its ability to perform MP2-level AIMD simulations on systems with thousands of atoms, leveraging novel algorithms optimized for GPU architectures (NVIDIA).

**Scientific domain**: High-performance quantum chemistry, Ab Initio Molecular Dynamics (AIMD)  
**Target user community**: HPC users, researchers needing large-scale accurate dynamics (MP2)

## Theoretical Methods
- Hartree-Fock (HF)
- Density Functional Theory (DFT)
- Second-order Møller-Plesset perturbation theory (MP2)
- Ab Initio Molecular Dynamics (AIMD)
- Resolution of Identity (RI) approximations
- GPU-accelerated algorithms

## Capabilities (CRITICAL)
- GPU-native implementation (CUDA)
- Extreme scalability (Summit, Frontier, Aurora scales)
- Large-scale MP2 calculations (1000+ atoms)
- Long-timescale AIMD at MP2 level
- Energy conservation in dynamics
- High floating-point efficiency
- Massively parallel execution

## Key Strengths

### GPU Optimization:
- Built from scratch for GPUs
- High percent peak flop utilization
- Minimal CPU-GPU transfer
- optimized Tensor contractions
- Scalable to thousands of GPUs

### MP2 Dynamics:
- Accurate electron correlation
- Dispersion inclusion via MP2
- Feasible for large biological systems
- Beyond DFT accuracy for dynamics

### Extreme Scale:
- Linear scaling or low-prefactor algorithms
- Handles 1000-2000 atoms at MP2 level
- Gordon Bell Prize performance
- State-of-the-art HPC

## Inputs & Outputs
- **Input formats**:
  - EXESS input format
  - PDB/XYZ coordinates
  - Basis set library inputs
- **Output data types**:
  - Energies (HF, MP2)
  - Forces
  - Trajectories (XYZ/DCD)
  - Restart files
  - Performance metrics

## Interfaces & Ecosystem
- **HPC Systems**: Designed for OLCF/ALCF supercomputers (Summit, Frontier, Aurora)
- **NVIDIA Integration**: Optimized for A100/H100 GPUs using CUDA and cuBLAS
- **Analysis Tools**: Standard trajectory analysis (VMD, MDAnalysis)
- **Input Generation**: Minimal input scripts, compatible with standard formats
- **Output Parsing**: Standard text output, easily parsable performance logs

## Advanced Features

### MP2-AIMD:
- On-the-fly forces
- Conserved energy dynamics
- Solvated systems
- Chemical reactions in solution

### Algorithms:
- Rank-reduced operations
- Mixed precision utilization
- Asynchronous task scheduling
- Distributed memory management

## Performance Characteristics
- **Speed**: Orders of magnitude faster than CPU codes for MP2
- **Accuracy**: MP2/CBS limit capabilities
- **System size**: 1000+ atoms (MP2)
- **Memory**: GPU memory constrained (managed)
- **Parallelization**: Multi-node Multi-GPU (MPI+CUDA)

## Computational Cost
- **MP2**: Conventionally O(N^5), EXESS optimized
- **Dynamics**: Feasible ps/ns scales
- **Hardware**: High-end GPU clusters required
- **Efficiency**: High FLOP/watt

## Limitations & Known Constraints
- **Availability**: Not open source (Waitlist)
- **Hardware**: Requires NVIDIA GPUs
- **Features**: Focused on energy/forces (MP2), less property analysis
- **Documentation**: Limited public docs

## Comparison with Other Codes
- **vs CP2K**: EXESS focuses on MP2, CP2K on DFT
- **vs GAMESS**: EXESS GPU-native, faster for large MP2
- **vs TeraChem**: Both GPU, EXESS targets HPC/MP2 scale
- **vs Psi4**: EXESS is HPC dynamics focused
- **Unique strength**: Large-scale MP2 dynamics on GPUs

## Application Areas

### Biochemistry:
- **Enzyme Reactions**: MP2 accuracy for reaction mechanisms in large enzymes
- **Solvation Dynamics**: Accurate description of solvation shells including dispersion
- **Ligand Binding**: Free energy calculations with correlated methods
- **Conformational Ensembles**: Sampling complex landscapes with high accuracy

### Materials Science:
- **Liquid Structures**: Reliable radial distribution functions from MP2
- **Interfacial Chemistry**: Solid-liquid interfaces with accurate electronic structure
- **Nanoparticles**: Dynamics of metallic and semiconductor clusters
- **Battery Electrolytes**: Solvation structures and transport mechanisms

## Best Practices
### System Setup:
- **Pre-equilibration**: thorough equilibration with classical MD before switching to EXESS
- **Basis Sets**: Use standard correlation-consistent basis sets (cc-pVDZ/TZ)
- **Geometry**: Ensure clean starting structures to avoid large initial forces

### Hardware Utilization:
- **GPU Selection**: Target A100 or H100 nodes for maximum efficiency
- **Memory Management**: Monitor GPU memory usage for large basis sets
- **Scaling**: Test scaling on small number of nodes before full production run

### Simulation Parameters:
- **Timestep**: Use appropriate timestep (0.5-1.0 fs) for AIMD
- **Thermostats**: Standard thermostats (Nose-Hoover) available
- **Restart frequency**: Write restarts frequently due to HPC time limits

## Community and Support
- Barca group (Australian National University)
- HPC centers (OLCF, etc.)
- Gordon Bell community
- Academic collaborations

## Verification & Sources
**Primary sources**:
1. Homepage: https://barcagrp.com/exess/
2. Gordon Bell Prize 2024 announcements
3. Barca group publications (JCTC, etc.)

**Confidence**: VERIFIED
- Status: Active HPC code
- Recognition: Gordon Bell Prize
- Existence: Confirmed via ANU/OLCF
