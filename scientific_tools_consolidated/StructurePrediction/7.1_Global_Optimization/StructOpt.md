# StructOpt

## Overview
StructOpt is a modular structure optimization suite designed for materials with complicated structures. It identifies atomic structures that are energetically stable and consistent with experimental data by combining multiple fitness criteria.

## Theoretical Basis
- Genetic algorithm optimization
- Multi-objective fitness evaluation
- Experimental data fitting (STEM, FEM)
- Energy minimization
- Structure fingerprinting

## Key Capabilities
- Multi-objective structure optimization
- Experimental data integration (STEM, FEM)
- Modular fitness functions
- Nanoparticle and defect optimization
- Amorphous structure prediction

**Sources**: Comp. Mater. Sci. 156, 204 (2019)

## Key Strengths

### Experimental Integration:
- STEM image fitting
- FEM data matching
- Multi-data fusion

### Modularity:
- Pluggable fitness functions
- Extensible architecture
- Custom objectives

### Applications:
- Nanoparticles
- Amorphous materials
- Defect structures

## Inputs & Outputs
- **Input formats**: Structure files, experimental data
- **Output data types**: Optimized structures, fitness scores

## Interfaces & Ecosystem
- **Calculators**: LAMMPS, VASP
- **Experimental**: STEM simulation, FEM
- **Analysis**: Structure comparison

## Workflow and Usage
1. Define system and constraints
2. Configure fitness functions
3. Provide experimental data (optional)
4. Run genetic algorithm
5. Analyze Pareto-optimal structures

## Performance Characteristics
- Depends on fitness evaluation cost
- Parallelizable populations
- Efficient for medium systems

## Computational Cost
- Calculator-dependent
- Experimental fitting adds overhead
- Parallelization helps

## Best Practices
- Use appropriate fitness weights
- Include experimental constraints
- Monitor population diversity
- Validate final structures

## Limitations & Known Constraints
- Complex setup for multi-objective
- Requires experimental data for full benefit
- Learning curve

## Application Areas
- Nanoparticle structure determination
- Amorphous-nanocrystal composites
- Defect structure optimization
- Experimental structure refinement

## Comparison with Other Codes
- **vs USPEX**: StructOpt experimental-focused
- **vs AGOX**: Different optimization targets
- **Unique strength**: Experimental data integration, multi-objective

## Community and Support
- Open-source (GitHub)
- Academic development (UW-Madison)
- Documentation available

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/uw-cmg/StructOpt_modular
2. Publication: Comp. Mater. Sci. 156, 204 (2019)
3. Documentation: https://structopt.readthedocs.io/

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE (GitHub)
- Documentation: AVAILABLE
- Source: OPEN
- Development: MAINTAINED
- Applications: Structure optimization with experimental data
