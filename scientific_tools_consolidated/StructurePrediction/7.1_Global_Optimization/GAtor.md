# GAtor

## Overview
GAtor is a massively parallel, first-principles genetic algorithm (GA) for molecular crystal structure prediction. Written in Python, it interfaces with the FHI-aims code for local optimizations and energy evaluations using dispersion-inclusive DFT.

## Theoretical Basis
- Genetic algorithm optimization
- First-principles DFT energy evaluation
- Dispersion-inclusive functionals (vdW-DF, TS, MBD)
- Evolutionary niching with machine learning clustering
- Crossover and mutation operators for molecular crystals

## Key Capabilities
- Molecular crystal structure prediction
- Massively parallel execution
- Machine learning-based niching
- Multiple fitness evaluation schemes
- Flexible breeding operators

**Sources**: J. Chem. Theory Comput. 14, 2246 (2018)

## Key Strengths

### Molecular Crystals:
- Designed for organic/molecular systems
- Handles conformational flexibility
- Proper treatment of dispersion

### Parallelization:
- Massively parallel execution
- Efficient HPC utilization
- Scales to large populations

### ML Integration:
- Clustering for evolutionary niching
- Diversity maintenance
- Efficient exploration

## Inputs & Outputs
- **Input formats**: Molecular geometry, GA parameters
- **Output data types**: Crystal structures, energies, population history

## Interfaces & Ecosystem
- **DFT**: FHI-aims (primary interface)
- **Dispersion**: vdW-DF, TS, MBD methods
- **Analysis**: Structure comparison, clustering

## Workflow and Usage
1. Prepare molecular geometry
2. Configure GA parameters
3. Set up FHI-aims interface
4. Run parallel GA search
5. Analyze converged structures

## Performance Characteristics
- Excellent parallel scaling
- Efficient for molecular crystals
- ML-assisted diversity maintenance

## Computational Cost
- DFT-limited (FHI-aims calculations)
- Parallelization reduces wall time
- Population size affects cost

## Best Practices
- Use appropriate dispersion correction
- Tune population size for system
- Enable niching for diversity
- Validate with experimental data

## Limitations & Known Constraints
- Requires FHI-aims license
- Computationally expensive (DFT)
- Focused on molecular crystals

## Application Areas
- Organic crystal structure prediction
- Pharmaceutical polymorph screening
- Molecular materials design
- Crystal engineering

## Comparison with Other Codes
- **vs USPEX**: GAtor molecular-focused, USPEX more general
- **vs Genarris**: GAtor includes optimization, Genarris generation only
- **vs MGAC**: Similar approach, different implementation
- **Unique strength**: First-principles GA, ML niching, molecular crystals

## Community and Support
- Academic development (CMU)
- Published methodology
- Available upon request

## Verification & Sources
**Primary sources**:
1. Publication: Curtis et al., J. Chem. Theory Comput. 14, 2246 (2018)
2. arXiv: https://arxiv.org/abs/1802.08602

**Secondary sources**:
1. FHI-aims documentation
2. Molecular crystal CSP literature

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Documentation: AVAILABLE (paper)
- Source: ACADEMIC (request)
- Development: ACTIVE
- Applications: Molecular crystal structure prediction
