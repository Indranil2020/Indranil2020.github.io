# BigDFT

## Official Resources
- Homepage: https://bigdft.org/
- Documentation: https://l_sim.gitlab.io/bigdft-suite/
- Source Repository: https://gitlab.com/l_sim/bigdft-suite
- PyBigDFT: https://pypi.org/project/BigDFT/
- License: GNU General Public License v2.0

## Overview
BigDFT is a DFT code using Daubechies wavelets as a basis set, providing systematic convergence and efficient treatment of systems in vacuum, surfaces, and periodic systems. It features linear-scaling capabilities, excellent support for massively parallel calculations, early GPU adoption, and a comprehensive Python interface (PyBigDFT) for workflow automation.

**Scientific domain**: Molecules, nanostructures, surfaces, linear-scaling calculations, biological systems  
**Target user community**: Researchers studying isolated systems, surfaces, and needing systematic basis convergence and Python-driven workflows

## Theoretical Methods
- Density Functional Theory (DFT)
- Daubechies wavelet basis sets
- Systematic basis set convergence
- Norm-conserving and HGH pseudopotentials
- LDA, GGA functionals (via LibXC)
- Hybrid functionals
- van der Waals corrections
- DFT+U for correlated systems
- Linear-scaling DFT (O(N) method)
- Time-Dependent DFT (in development)
- Poisson solver for arbitrary boundary conditions

## Capabilities (CRITICAL)
- Ground-state electronic structure
- Systematic convergence with single parameter
- Isolated molecules (no spurious interactions)
- Surfaces and low-dimensional systems
- Periodic systems (1D, 2D, 3D)
- Linear-scaling DFT for large systems (thousands of atoms)
- Geometry optimization and transition states
- Molecular dynamics (NVE, NVT, NPT)
- Band structure and DOS
- Forces and stress tensors
- Massively parallel calculations (thousands of cores)
- GPU acceleration (CUDA, OpenCL)
- Adaptive grid refinement
- Fragment-based calculations
- Implicit electrostatic solvents
- Fragmentation analysis (charges, dipoles, interactions)

**Sources**: Official BigDFT documentation, cited in 6/7 source lists

## Key Strengths

### Daubechies Wavelet Basis:
- Systematic convergence (single cutoff parameter)
- Compact support (localized)
- Adaptive mesh refinement
- Mathematical rigor
- No basis set superposition error

### Flexible Boundary Conditions:
- Free boundary (isolated molecules)
- Surface boundary (slabs)
- Wire boundary (1D periodic)
- Periodic (3D crystals)
- Mixed boundary conditions

### Linear-Scaling Approach:
- Fragment-based O(N) method
- Thousands of atoms feasible
- Automatic system partitioning
- Localized description from wavelets
- Production large-scale calculations

### GPU Acceleration:
- Early CUDA adoption
- OpenCL support
- Hybrid CPU/GPU execution
- Significant speedups
- Modern HPC ready

### PyBigDFT Python Interface:
- High-level Python API
- Workflow automation
- Jupyter notebook support
- Post-processing and analysis
- Interoperability with other tools

## Inputs & Outputs
- **Input formats**:
  - Input files (YAML format)
  - XYZ coordinate files
  - Pseudopotential files
  - Python API for input generation
  
- **Output data types**:
  - YAML output files
  - Energies and forces
  - Optimized structures
  - Wavefunction data
  - Density files
  - DOS outputs
  - Fragment analysis

## Interfaces & Ecosystem
- **Python integration**:
  - PyBigDFT - comprehensive Python package
  - Jupyter notebook support
  - High-level workflow API
  - Post-processing tools
  - Scriptable automation
  
- **Framework integrations**:
  - ASE interface
  - Babel
  - RDKit
  - XTB
  - OpenMM
  - PSI4
  - DFTB
  - MRChem
  
- **Visualization**:
  - v_sim - visualization tool
  - Compatible with standard tools
  - Density visualization
  - Fragment visualization
  
- **Linear-scaling**:
  - Fragment approach
  - Excellent parallel scaling
  - Automatic fragmentation

## Advanced Features

### Fragmentation Analysis:
- Automatic system partitioning
- Fragment charges and dipoles
- Inter-fragment interactions
- Energy decomposition
- Locality analysis

### BioQM Module:
- Biological system analysis
- Enzyme studies
- Protein-ligand interactions
- Specialized biomolecular tools

### Implicit Solvation:
- Electrostatic solvation models
- Environment effects
- Solvation free energies
- Boundary condition handling

### Workflow Automation:
- Python-driven calculations
- High-throughput capable
- Result parsing and analysis
- Plotting and visualization
- Custom analysis scripts

### Interoperability:
- ASE atoms and calculators
- Babel molecular formats
- RDKit chemistry tools
- XTB semi-empirical
- OpenMM molecular dynamics

## Performance Characteristics
- **Speed**: Efficient wavelet implementation
- **Accuracy**: Systematic with single parameter
- **System size**: Thousands of atoms with O(N)
- **Memory**: Efficient for wavelets
- **Parallelization**: Excellent GPU and MPI scaling

## Computational Cost
- **DFT**: Efficient wavelet evaluation
- **O(N)**: Linear scaling achieved
- **GPU**: Significant acceleration
- **Large systems**: Fragment approach efficient
- **Typical**: Competitive with major codes

## Limitations & Known Constraints
- **Wavelet basis**: Less familiar than plane-waves or orbitals
- **Pseudopotentials**: Requires specific HGH format
- **Hybrid functionals**: Limited implementation
- **Community**: Smaller than major codes
- **Documentation**: Good but evolving
- **Learning curve**: Wavelet methods require understanding
- **k-point sampling**: Best for systems where Γ-point sufficient
- **Installation**: Requires compilation and libraries
- **Platform**: Primarily Linux/Unix

## Comparison with Other Codes
- **vs VASP/QE**: BigDFT wavelets, plane-wave codes pseudopotentials
- **vs CONQUEST/ONETEP**: All O(N), different basis technologies
- **vs NWChem**: BigDFT more focused on wavelets and large systems
- **vs Gaussian**: BigDFT materials focus, Gaussian molecular
- **Unique strength**: Daubechies wavelets, systematic convergence, GPU acceleration, PyBigDFT Python interface, fragmentation analysis

## Application Areas

### Biological Systems:
- Enzymes and proteins
- Drug-target interactions
- QM studies of active sites
- Large biomolecular complexes
- BioQM specialized analysis

### Nanomaterials:
- Nanoparticles
- Clusters
- Functionalized surfaces
- Quantum dots
- Carbon nanostructures

### Isolated Molecules:
- No periodic image interactions
- Accurate for finite systems
- Molecular properties
- Reaction barriers
- Conformational analysis

### Surfaces and Interfaces:
- Surface chemistry
- Adsorption studies
- Interface electronic structure
- Electrostatic boundary handling

### OLED Materials:
- Organic electronics
- Excited state properties
- Large molecular systems
- Materials design

## Best Practices

### Grid Convergence:
- Single hgrid parameter
- Systematic convergence
- Test total energy convergence
- Balance accuracy vs cost

### Boundary Conditions:
- Free for molecules
- Surface for slabs
- Periodic for crystals
- Match to physical system

### Fragment Calculations:
- Appropriate fragment size
- Monitor interaction energies
- Validate fragmentation
- Use for very large systems

### Python Workflow:
- Use PyBigDFT for complex workflows
- Jupyter for interactive analysis
- Script high-throughput
- Leverage interoperability

## Community and Support
- Open-source GPL v2
- Active GitLab development (20,000+ commits)
- Documentation and tutorials
- PyBigDFT on PyPI
- Video tutorials available
- Regular releases (13 on GitLab)

## Verification & Sources
**Primary sources**:
1. Official website: https://bigdft.org/
2. Documentation: https://l_sim.gitlab.io/bigdft-suite/
3. GitLab repository: https://gitlab.com/l_sim/bigdft-suite
4. L. Genovese et al., J. Chem. Phys. 129, 014109 (2008) - Daubechies wavelets
5. S. Mohr et al., J. Chem. Phys. 140, 204110 (2014) - BigDFT linear-scaling

**Secondary sources**:
1. BigDFT tutorials and workshops
2. PyBigDFT documentation
3. Published applications
4. Confirmed in 6/7 source lists (claude, g, gr, k, m, q)

**Confidence**: CONFIRMED - Appears in 6 of 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: ACCESSIBLE
- Source code: OPEN (GitLab, GPL v2)
- Community support: Active (developers, documentation)
- Academic citations: >300 (main papers)
- Active development: Regular releases (13 on GitLab), 20,000+ commits
- Specialized strength: Daubechies wavelets, systematic convergence, GPU acceleration, PyBigDFT Python ecosystem, fragmentation analysis, BioQM
