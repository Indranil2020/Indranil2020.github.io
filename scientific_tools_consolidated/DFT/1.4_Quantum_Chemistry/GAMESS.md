# GAMESS (US)

## Official Resources
- Homepage: https://www.msg.chem.iastate.edu/gamess/
- Documentation: https://www.msg.chem.iastate.edu/gamess/documentation.html
- Source Repository: Available with license
- License: Free for academic and commercial use (registration required)

## Overview
GAMESS (General Atomic and Molecular Electronic Structure System) is a comprehensive ab initio quantum chemistry package developed at Iowa State University. Free for all users, GAMESS provides extensive capabilities from Hartree-Fock to high-level correlated methods, with particular strengths in excited states, solvation, and QM/MM calculations. It is one of the most widely distributed quantum chemistry codes worldwide.

**Scientific domain**: Quantum chemistry, molecular modeling, excited states, solvation  
**Target user community**: Computational chemists across academia and industry worldwide

## Theoretical Methods
- Hartree-Fock (RHF, UHF, ROHF)
- Density Functional Theory (DFT)
- LDA, GGA, meta-GGA, hybrid functionals
- MP2, MP3, MP4
- Coupled cluster (CCSD, CCSD(T), CR-CC)
- Multi-reference methods (MCSCF, MRMP2, MRCI)
- Complete active space (CASSCF, CASPT2)
- Time-Dependent DFT (TDDFT)
- Configuration interaction (CI, CISD)
- Solvation models (PCM, SMD, EFP)
- QM/MM methods
- Effective Fragment Potential (EFP)
- Relativistic methods (DKH, RESC)
- Spin-orbit coupling
- Fragment molecular orbital (FMO)

## Capabilities (CRITICAL)
- Ground-state electronic structure
- Geometry optimization and transition states
- Intrinsic reaction coordinate (IRC)
- Vibrational frequencies and thermochemistry
- Excited states (TDDFT, MCSCF, CI, EOM-CC)
- Conical intersections
- Surface hopping dynamics
- NMR chemical shifts and spin-spin coupling
- Optical rotation and circular dichroism
- Vibrational circular dichroism (VCD)
- Solvation free energies
- QM/MM calculations
- Fragment molecular orbital (large systems)
- Effective Fragment Potential
- Molecular dynamics (Born-Oppenheimer)
- Analytic gradients for many methods
- Parallel execution (MPI, OpenMP)
- GPU acceleration (limited)

**Sources**: Official GAMESS documentation (https://www.msg.chem.iastate.edu/gamess/), confirmed in multiple source lists

## Key Strengths

### Free Availability:
- No cost for anyone
- Academic and commercial use
- Worldwide distribution
- Registration only
- Source code access

### Comprehensive Methods:
- HF to high-level correlation
- Multi-reference capabilities
- Excited state methods
- Solvation models
- QM/MM integration

### Effective Fragment Potential:
- Unique EFP method
- Fast solvent treatment
- Explicit solvent without full QM
- Polarization effects
- Accurate interactions

### Fragment Molecular Orbital:
- Linear-scaling for large systems
- Protein calculations
- Thousands of atoms
- Accurate energies
- Parallelizable

### Community:
- Large user base
- Long development history
- Well-tested
- Extensive documentation
- Active forum

## Inputs & Outputs
- **Input formats**:
  - Text-based input file
  - $CONTRL, $BASIS, etc. groups
  - Z-matrix or Cartesian coordinates
  - Internal coordinate definitions
  
- **Output data types**:
  - Detailed text output
  - Energies and gradients
  - Molecular orbitals
  - Excited state information
  - Vibrational modes
  - Trajectory files

## Interfaces & Ecosystem
- **GUIs**:
  - MacMolPlt (molecular builder/viewer)
  - wxMacMolPlt
  - Third-party interfaces (Avogadro, ChemCraft)
  
- **Visualization**:
  - MacMolPlt for orbitals
  - Standard visualization tools
  - Trajectory analysis
  
- **Job Management**:
  - GamessQ (GUI for job submission)
  - Batch scripts
  - Cluster integration
  
- **Parallelization**:
  - MPI for distributed memory
  - DDI (Distributed Data Interface)
  - OpenMP for shared memory
  - Hybrid parallelization

## Workflow and Usage

### Input File Structure:

```
 $CONTRL SCFTYP=RHF RUNTYP=OPTIMIZE $END
 $SYSTEM TIMLIM=525600 MEMORY=1000000 $END
 $BASIS  GBASIS=N31 NGAUSS=6 NDFUNC=1 $END
 $STATPT OPTTOL=0.0001 NSTEP=50 $END
 $DATA
Water molecule
C1
O  8.0  0.0  0.0  0.0
H  1.0  0.0  0.0  1.0
H  1.0  0.0  1.0  0.0
 $END
```

### Running GAMESS:
```bash
rungms water 01 1
# version 01, 1 processor

rungms water 01 4
# 4 processors
```

### Common Calculations:
- Energy: RUNTYP=ENERGY
- Optimization: RUNTYP=OPTIMIZE
- Frequencies: RUNTYP=HESSIAN
- IRC: RUNTYP=IRC
- TDDFT: TDDFT=EXCITE

## Advanced Features

### Effective Fragment Potential:
- Fast solvent treatment
- Polarizable force field
- Coulomb, polarization, dispersion
- Water, organic solvents
- Biological systems

### Fragment Molecular Orbital:
- Linear-scaling method
- Fragment decomposition
- Large biomolecules
- Protein-ligand binding
- Massively parallel

### Multi-Reference Methods:
- MCSCF for multi-configurational
- CASSCF for active space
- MRMP2 for dynamic correlation
- MRCI for high accuracy
- Conical intersections

### Surface Hopping:
- Non-adiabatic dynamics
- Excited state trajectories
- Fewest switches algorithm
- Photochemistry
- UV spectroscopy

### Solvation:
- PCM (Polarizable Continuum Model)
- SMD solvation model
- EFP explicit solvent
- COSMO
- Multiple approaches

## Performance Characteristics
- **Speed**: Competitive
- **Scaling**: Good with MPI
- **Memory**: Moderate to high
- **System size**: Up to ~500 atoms (standard); thousands (FMO)
- **Parallelization**: Efficient MPI implementation

## Computational Cost
- **HF/DFT**: Standard scaling
- **MP2**: Manageable for medium systems
- **CCSD(T)**: Expensive, high accuracy
- **MCSCF**: Expensive, small active spaces
- **FMO**: Efficient for very large systems
- **EFP**: Very fast for solvents

## Limitations & Known Constraints
- **Learning curve**: Steep input format
- **GUI**: Limited official GUI support
- **Modern features**: Fewer than commercial codes
- **Documentation**: Comprehensive but dense
- **Parallelization**: Good but not cutting-edge
- **GPU**: Limited GPU support
- **Platform**: Linux, macOS, Windows

## Comparison with Other Codes
- **vs Gaussian**: GAMESS free, Gaussian commercial; similar capabilities
- **vs NWChem**: Both free, different implementations
- **vs ORCA**: ORCA more modern, GAMESS more established
- **vs Psi4**: Psi4 more modern, GAMESS broader methods
- **Unique strength**: Free, EFP method, FMO for large systems, comprehensive multi-reference

## Application Areas

### Molecular Chemistry:
- Organic molecules
- Reaction mechanisms
- Thermochemistry
- Conformational analysis
- Structure determination

### Excited States:
- UV-Vis spectroscopy
- Fluorescence/phosphorescence
- Photochemistry
- Charge transfer
- Non-adiabatic dynamics

### Biochemistry:
- Protein-ligand binding (FMO)
- Enzyme mechanisms
- Drug design
- Biomolecular properties
- Large systems

### Solvation:
- Solution-phase chemistry
- Solvent effects
- Free energies
- Explicit solvent (EFP)
- Implicit models

## Best Practices

### Input Preparation:
- Use MacMolPlt for building
- Check input syntax
- Start with simple calculations
- Test basis set convergence
- Verify SCF convergence

### Method Selection:
- HF/DFT for large systems
- MP2 for moderate correlation
- CCSD(T) for benchmarks
- MCSCF for multi-reference
- FMO for very large systems

### Basis Sets:
- 6-31G(d) for quick tests
- 6-311G(d,p) for publication
- cc-pVTZ for high accuracy
- Augmented for anions/excited states

### Convergence:
- Appropriate SCF settings
- Good initial geometry
- Symmetry when applicable
- Check for convergence issues
- Use stability analysis

### Parallelization:
- Test scaling efficiency
- Balance processors/memory
- Use MPI for distributed
- Hybrid for large nodes

## Community and Support
- Free worldwide distribution
- Active user forum
- Mailing lists
- Comprehensive manual
- Regular updates
- Iowa State development

## Educational Resources
- Detailed manual (>500 pages)
- Example inputs
- Tutorial workshops
- Published papers
- User forum knowledge
- Video tutorials (community)

## Development
- Iowa State University
- Mark Gordon research group
- Active development since 1980s
- Regular releases
- Community contributions
- Open collaboration

## Historical Significance
- One of oldest quantum chemistry codes
- Pioneered free distribution
- Trained generations of chemists
- Extensive method development
- Worldwide impact

## Special Features

### EFP Method:
- Unique to GAMESS
- Fast, accurate solvent
- Polarizable model
- Library of fragments
- QM/EFP hybrid

### FMO Method:
- Large biomolecule capability
- Linear scaling
- Fragment analysis
- Binding energies
- Interaction energies

### DDI:
- Distributed Data Interface
- Efficient parallelization
- Shared memory emulation
- Scalable

## Verification & Sources
**Primary sources**:
1. Official website: https://www.msg.chem.iastate.edu/gamess/
2. Documentation: https://www.msg.chem.iastate.edu/gamess/documentation.html
3. M. W. Schmidt et al., J. Comput. Chem. 14, 1347 (1993) - GAMESS overview
4. M. S. Gordon and M. W. Schmidt, Adv. Electron. Struct. Theory: GAMESS (2005)

**Secondary sources**:
1. GAMESS manual and documentation
2. Published studies using GAMESS (>20,000 citations)
3. User forum discussions
4. Confirmed in multiple source lists

**Confidence**: VERIFIED - Well-established, widely used code

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE (detailed manual)
- Software: Free (registration required)
- Community support: Active forum, mailing lists
- Academic citations: >25,000
- Active development: Regular releases from Iowa State
- Specialized strength: Free availability, EFP method, FMO for large systems, comprehensive methods, multi-reference capabilities, worldwide distribution
