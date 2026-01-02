# USPEX

## Official Resources
- Homepage: https://uspex-team.org/
- Documentation: https://uspex-team.org/en/uspex/documentation
- Source Repository: Available to registered users
- License: Free for academic use (registration required)

## Overview
USPEX (Universal Structure Predictor: Evolutionary Xtallography) is a method and code for crystal structure prediction using evolutionary algorithms. Developed by Artem Oganov and collaborators, it is the leading tool for predicting stable crystal structures, interfaces, surfaces, and nanoparticles from scratch, requiring only chemical composition. It has led to numerous experimental discoveries of new materials.

**Scientific domain**: Crystal structure prediction, materials discovery, evolutionary algorithms  
**Target user community**: Materials scientists, crystallographers seeking novel structures and phases

## Theoretical Methods
- Evolutionary algorithm for structure prediction
- Variation operators (heredity, mutation, permutation)
- Fitness function based on enthalpy/free energy
- Local optimization via external codes
- Particle swarm optimization
- Metadynamics for free energy landscapes
- Fingerprint-based structure analysis
- Space group symmetry constraints
- Variable composition searches
- Interface and surface structure prediction
- Fingerprint-based structure comparison
- Fitness function based on enthalpy/free energy from DFT
- Multi-objective optimization

## Capabilities (CRITICAL)
- Crystal structure prediction for bulk materials (0D to 3D periodicity)
- Variable composition searches
- Surface and interface structure prediction
- Nanoparticle and cluster structure prediction
- 2D materials structure prediction
- Molecular crystal structure prediction
- Property-targeted structure search
- Hardness maximization
- Interface to multiple DFT codes (VASP, CASTEP, GULP, Quantum ESPRESSO, CP2K, SIESTA, ABINIT, CRYSTAL)
- Interface to classical potentials (GULP, LAMMPS)
- Interface to machine learning potentials
- Parallel execution with distributed calculations
- Constraint-based searches (symmetry, coordination)
- Anti-seeding to avoid known structures

**Sources**: Official USPEX website, publications, cited in 7/7 source lists

## Inputs & Outputs
- **Input formats**:
  - INPUT.txt (main USPEX input file)
  - INCAR_1, INCAR_2, etc. (for VASP interface)
  - Seed structures (optional initial guesses)
  - Constraint files for specific searches
  
- **Output data types**:
  - gatheredPOSCARS (all generated structures)
  - BESTgatheredPOSCARS (best structures per generation)
  - Individuals files (structure database)
  - results.txt (evolution progress)
  - Convex hull data
  - Hardness predictions
  - Phonon stability checks

## Interfaces & Ecosystem
- **DFT code interfaces** (verified):
  - VASP - primary interface, most tested
  - CASTEP - supported
  - Quantum ESPRESSO - supported
  - CP2K - supported
  - SIESTA - supported
  - ABINIT - supported
  - CRYSTAL - supported
  - GULP - for classical potentials
  - LAMMPS - for MD and classical potentials
  
- **Post-processing tools**:
  - Built-in convex hull analysis
  - Structure visualization utilities
  - Phonon stability checking
  - Hardness prediction modules
  
- **Machine learning integration**:
  - Can interface with ML potentials for acceleration
  - GAP, moment tensor potentials supported

## Limitations & Known Constraints
- **Computational cost**: Structure prediction extremely expensive; requires hundreds to thousands of DFT calculations
- **Registration required**: Not fully open-source; requires registration for download
- **DFT dependency**: Quality depends on underlying DFT calculations; expensive for hybrid functionals
- **System size**: Limited to ~20-50 atoms per unit cell for ab initio; larger with classical potentials
- **Convergence**: No guarantee of finding global minimum; multiple runs often needed
- **Initial population size**: Requires careful tuning for specific systems
- **Success rate**: Higher for simple systems; complex multi-component systems more challenging
- **Learning curve**: Requires understanding of evolutionary algorithms and DFT
- **Parallelization**: Requires job submission system for distributed DFT calculations
- **Post-processing**: Manual analysis of large numbers of structures can be tedious

## Verification & Sources
**Primary sources**:
1. Official website: https://uspex-team.org/
2. C. W. Glass et al., Comput. Phys. Commun. 175, 713 (2006) - Original USPEX paper
3. A. R. Oganov and C. W. Glass, J. Chem. Phys. 124, 244704 (2006) - Evolutionary algorithm
4. A. O. Lyakhov et al., Comput. Phys. Commun. 184, 1172 (2013) - Variable composition searches
5. Q. Zhu et al., Phys. Rev. Lett. 106, 145501 (2011) - Surfaces and nanoparticles

**Secondary sources**:
1. USPEX tutorials and user manual
2. Workshop materials and examples
3. Published structure predictions using USPEX (numerous papers)
4. Confirmed in 7/7 source lists (claude, g, gr, k, m, q, z)

**Confidence**: CONFIRMED - Appears in all 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: ACCESSIBLE (requires registration)
- User community: Large and active
- Academic citations: >1,500 (main papers)
- Success stories: Multiple novel materials predicted and synthesized
- Active development: Regular updates and new features
