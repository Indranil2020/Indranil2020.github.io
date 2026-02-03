# AGOX

## Overview
AGOX (Atomistic Global Optimization X) is a Python package for global optimization of atomistic structures. It interfaces with ASE and supports any ASE-compatible calculator, making it highly flexible for various optimization tasks.

## Theoretical Basis
- Basin hopping algorithm
- Genetic algorithms
- Particle swarm optimization
- Machine learning surrogate models
- Bayesian optimization
- Modular sampler/evaluator architecture

## Key Capabilities
- Multiple global optimization algorithms
- ASE calculator compatibility
- Machine learning acceleration
- Modular and extensible design
- Cluster, surface, and bulk optimization

**Sources**: arXiv:2204.01451, AGOX documentation

## Key Strengths

### Flexibility:
- Any ASE calculator
- Modular architecture
- Easy customization

### Algorithms:
- Basin hopping, GA, PSO
- Bayesian optimization
- ML surrogate models

### Applications:
- Clusters, surfaces, bulk
- Defects, interfaces
- Nanoparticles

## Inputs & Outputs
- **Input formats**: ASE Atoms objects, configuration files
- **Output data types**: Optimized structures, trajectories, databases

## Interfaces & Ecosystem
- **Calculators**: Any ASE-compatible (VASP, GPAW, EMT, ML potentials)
- **ASE**: Full integration
- **Databases**: ASE database support

## Workflow and Usage
1. Define system (composition, constraints)
2. Select calculator and algorithm
3. Configure AGOX settings
4. Run optimization
5. Analyze results

## Performance Characteristics
- Depends on calculator choice
- ML acceleration available
- Efficient for clusters and surfaces

## Computational Cost
- Calculator-limited
- ML surrogate reduces evaluations
- Parallelizable

## Best Practices
- Use ML surrogate for expensive calculators
- Choose algorithm based on system
- Enable structure comparison
- Monitor convergence

## Limitations & Known Constraints
- Requires ASE knowledge
- Less specialized than dedicated CSP codes
- Documentation improving

## Application Areas
- Cluster structure optimization
- Surface reconstruction
- Nanoparticle structure
- Defect configurations
- Interface structures

## Comparison with Other Codes
- **vs ASE-GA**: AGOX more algorithms, more modular
- **vs GMIN**: AGOX Python/ASE-based, GMIN Fortran
- **vs CrySPY**: Different focus (clusters vs crystals)
- **Unique strength**: ASE integration, modular design, ML acceleration

## Community and Support
- Open-source (GPLv3)
- GitLab/GitHub repository
- Active development
- Documentation available

## Verification & Sources
**Primary sources**:
1. GitLab: https://agox.gitlab.io/agox/
2. GitHub: https://github.com/kimrojas/agox
3. arXiv: https://arxiv.org/abs/2204.01451

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN (GPLv3)
- Development: ACTIVE
- Applications: Global optimization, clusters, surfaces
