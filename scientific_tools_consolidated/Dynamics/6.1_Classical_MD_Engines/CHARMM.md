# CHARMM (Chemistry at HARvard Macromolecular Mechanics)

## Official Resources
- Homepage: https://www.charmm.org/
- Documentation: https://www.charmm.org/charmm/documentation/
- Source Repository: https://www.charmm.org/ (Distributed via license)
- License: Proprietary / Academic License

## Overview
CHARMM is a highly versatile and widely used molecular simulation program with broad application to many-particle systems. It has been developed for over three decades, primarily at Harvard University. It provides a vast set of tools for molecular mechanics, molecular dynamics, free energy calculations, and analysis of biomolecules, polymers, and liquids.

**Scientific domain**: Biomolecular simulation, molecular mechanics, drug design
**Target user community**: Biophysicists, computational chemists, structural biologists

## Theoretical Methods
- Classical Molecular Dynamics (NVE, NVT, NPT)
- Molecular Mechanics / Energy Minimization
- Free Energy Perturbation (FEP) and Thermodynamic Integration (TI)
- Implicit Solvent Models (GB, PB, EEF1)
- Normal Mode Analysis (NMA)
- QM/MM (Quantum Mechanics / Molecular Mechanics) interfaces
- Replica Exchange (REMD)
- Path sampling methods

## Capabilities (CRITICAL)
- Extensive force field support (CHARMM force fields for proteins, lipids, nucleic acids, carbohydrates, small molecules)
- DOMDEC (Domain Decomposition) for high-performance parallel execution
- OpenMM integration for GPU acceleration
- Advanced vibrational analysis (quasi-harmonic analysis)
- Correlation function analysis
- Monte Carlo modules
- Scriptable command interface

**Sources**: CHARMM website, J. Comput. Chem. 30, 1545 (2009)

## Key Strengths

### Force Fields:
- CHARMM force fields (C36, C36m)
- CGenFF for small molecules
- Drude polarizable
- Extensively validated

### Methods:
- Advanced free energy (FEP, TI)
- QM/MM capabilities
- Normal mode analysis
- Path sampling

### CHARMM-GUI:
- Powerful web interface
- System setup automation
- Membrane builder
- Ligand parameterization

## Inputs & Outputs
- **Input formats**: Input scripts (.inp), Topology files (.rtf), Parameter files (.prm), Coordinate files (.crd/.pdb)
- **Output data types**: Trajectories (.dcd), Output logs (.out), Restart files (.res)

## Interfaces & Ecosystem
- **CHARMM-GUI**: Powerful web-based graphical interface for system setup
- **VMD**: Visualization
- **OpenMM**: GPU acceleration interface
- **MMTSB Toolset**: Perl-based toolkit for enhanced sampling
- **BioExcel**: Integration in workflows

## Workflow and Usage
1. **System Setup**: Use CHARMM-GUI or manual scripts to generate PSF (protein structure file) and CRD files.
2. **Minimization**: Run energy minimization to relax steric clashes.
3. **Equilibration**: Heat the system and equilibrate density/pressure.
4. **Production**: Run long timescale MD.
5. **Analysis**: Use CHARMM analysis facilities or external tools (VMD, cpptraj).

## Performance Characteristics
- Highly optimized for CPU clusters using DOMDEC
- GPU acceleration available via OpenMM or BLaDE (Basic Lambda Dynamics Engine)
- Scaling depends on system size and parallelization scheme

## Computational Cost
- Good parallel scaling (DOMDEC)
- GPU via OpenMM/BLaDE
- Efficient for medium systems
- Overall: Good for method development

## Best Practices
- Use CHARMM-GUI for setup
- Validate force field parameters
- Use appropriate ensemble
- Check energy conservation
- Use DOMDEC for parallel runs

## Limitations & Known Constraints
- Commercial/academic license
- Steeper learning curve
- Less GPU-optimized than AMBER
- Complex scripting language

## Application Areas
- Protein folding and stability
- Ligand binding affinity
- Lipid membrane dynamics
- Nucleic acid interactions
- Enzyme catalysis (QM/MM)

## Comparison with Other Codes
- **vs AMBER**: CHARMM more methods, AMBER better GPU
- **vs GROMACS**: CHARMM more flexible, GROMACS faster
- **vs NAMD**: CHARMM more features, NAMD better scaling
- **Unique strength**: CHARMM force fields, CHARMM-GUI, QM/MM, method development

## Community and Support
- Managed by the CHARMM Development Project
- Academic and commercial licensing available
- Active forums and user community
- Annual workshops

## Verification & Sources
**Primary sources**:
1. Homepage: https://www.charmm.org/
2. Publication: Brooks et al., J. Comput. Chem. 30, 1545 (2009)

**Secondary sources**:
1. CHARMM-GUI tutorials
2. CHARMM force field publications
3. Extensive published applications (>20,000 citations)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: LICENSED
- Development: ACTIVE (Harvard/community)
- Applications: Standard biomolecular MD, force field development
