# PLUMED

## Official Resources
- Homepage: https://www.plumed.org/
- Documentation: https://www.plumed.org/doc
- Source Repository: https://github.com/plumed/plumed2
- License: GNU Lesser General Public License v3.0

## Overview
PLUMED is an open-source library for free energy calculations in molecular systems which works together with some of the most popular molecular dynamics engines. It performs enhanced sampling calculations (Metadynamics, Umbrella Sampling, etc.) and analyzes trajectories using a wide range of collective variables. PLUMED can be used as a plugin or a standalone analysis tool.

**Scientific domain**: Enhanced sampling, free energy calculations, molecular dynamics analysis  
**Target user community**: Computational chemists, biophysicists, materials scientists

## Theoretical Methods
- Metadynamics (standard, well-tempered, bias-exchange)
- Umbrella Sampling
- Steered Molecular Dynamics
- Replica Exchange (with bias)
- Variationally Enhanced Sampling
- Thermodynamic Integration
- Reweighting techniques

## Capabilities (CRITICAL)
- Enhanced sampling using collective variables (CVs)
- Huge library of CVs (distances, angles, RMSD, coordination numbers, etc.)
- On-the-fly bias potentials
- Post-processing and analysis of trajectories
- Interface with MD codes (GROMACS, LAMMPS, NAMD, AMBER, CP2K, QE, VASP, etc.)
- Multiple walker metadynamics
- Machine learning collective variables (via PyTorch/TensorFlow)

**Sources**: PLUMED documentation, Comp. Phys. Comm. 185, 608 (2014)

## Key Strengths

### Universality:
- Works with 20+ MD codes
- Plugin architecture
- Standalone analysis
- Consistent interface

### CV Library:
- 100+ collective variables
- Custom CVs via Python
- Machine learning CVs
- Extensible

### Methods:
- Metadynamics variants
- Umbrella sampling
- Steered MD
- Reweighting

## Inputs & Outputs
- **Input formats**: PLUMED input script (text), MD trajectories (pdb, xtc, trr, dcd)
- **Output data types**: COLVAR files (time series of CVs), HILLS files (bias potentials), grid files

## Interfaces & Ecosystem
- **MD Engines**: Works as a patch or plugin for GROMACS, LAMMPS, NAMD, AMBER, CP2K, Quantum ESPRESSO, VASP, i-PI, OpenMM, etc.
- **Python**: wrappers for analysis
- **VMD**: Plugin for visualization

## Workflow and Usage
1. Patch MD code with PLUMED (if not built-in)
2. Create `plumed.dat` defining CVs and bias
3. Run MD code with PLUMED flag: `gmx mdrun -plumed plumed.dat`
4. Analyze output: `plumed sum_hills --hills HILLS`

## Performance Characteristics
- Minimal overhead for standard CVs
- Scalable with MPI (depends on MD code integration)
- Efficient grid-based bias evaluation

## Computational Cost
- Minimal overhead for simple CVs
- Complex CVs add cost
- Efficient grid-based bias
- Overall: Low overhead

## Best Practices
- Start with simple CVs
- Validate CV choice with unbiased runs
- Check metadynamics convergence
- Use well-tempered metadynamics
- Reweight for unbiased estimates

## Limitations & Known Constraints
- Requires patching some MD codes
- CV choice is critical
- Convergence can be slow
- Learning curve for advanced methods

## Application Areas
- Protein folding and conformational changes
- Chemical reactions in solution
- Crystal nucleation and growth
- Phase transitions
- Drug binding affinity

## Comparison with Other Codes
- **vs Colvars**: PLUMED more methods, Colvars native to NAMD
- **vs SSAGES**: PLUMED plugin-based, SSAGES standalone
- **vs OpenPathSampling**: PLUMED bias-based, OPS path-based
- **Unique strength**: Universal plugin, huge CV library, active consortium

## Community and Support
- Open-source (LGPL v3)
- Very active mailing list
- Masterclass tutorials
- Annual meetings
- Consortium-led development

## Verification & Sources
**Primary sources**:
1. Homepage: https://www.plumed.org/
2. GitHub: https://github.com/plumed/plumed2
3. Publication: Comp. Phys. Comm. 185, 608 (2014); Nat. Methods 16, 670 (2019)

**Secondary sources**:
1. PLUMED Masterclass tutorials
2. PLUMED-NEST repository
3. Extensive published applications (>5000 citations)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN (GitHub)
- Development: ACTIVE (PLUMED Consortium)
- Applications: Enhanced sampling, metadynamics, free energy, MD plugin
