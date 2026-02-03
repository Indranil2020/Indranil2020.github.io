# ParetoCSP

## Overview
ParetoCSP is a crystal structure prediction algorithm that combines a multi-objective genetic algorithm (MOGA) with neural network interatomic potentials (M3GNet) to find energetically optimal crystal structures.

## Theoretical Basis
- Age-fitness Pareto genetic algorithm
- Multi-objective optimization
- M3GNet neural network potential
- Pareto front evolution
- Structure diversity maintenance

## Key Capabilities
- Multi-objective crystal structure prediction
- M3GNet integration for fast evaluation
- Age-fitness selection for diversity
- Pareto-optimal structure discovery
- Efficient global search

**Sources**: J. Mater. Inform. 4, 2 (2024)

## Key Strengths

### Multi-Objective:
- Age-fitness Pareto selection
- Diversity maintenance
- Multiple criteria optimization

### ML Acceleration:
- M3GNet potential
- Fast energy evaluation
- DFT-level accuracy

### Efficiency:
- Strong global search
- Reduced DFT calculations
- Competitive performance

## Inputs & Outputs
- **Input formats**: Chemical composition, constraints
- **Output data types**: Pareto-optimal structures, energies

## Interfaces & Ecosystem
- **ML Potentials**: M3GNet (primary)
- **Structure tools**: Pymatgen
- **Validation**: DFT codes

## Workflow and Usage
1. Define chemical composition
2. Configure GA parameters
3. Run multi-objective optimization
4. Extract Pareto front structures
5. Validate with DFT

## Performance Characteristics
- Fast with M3GNet
- Excellent global search
- Maintains structural diversity

## Computational Cost
- M3GNet evaluation: fast
- GA overhead: minimal
- DFT validation: optional

## Best Practices
- Use appropriate population size
- Enable age-fitness selection
- Validate top structures with DFT
- Check structural diversity

## Limitations & Known Constraints
- M3GNet accuracy limitations
- Requires validation for novel chemistries
- GA parameter tuning needed

## Application Areas
- Inorganic crystal structure prediction
- Materials discovery
- Phase prediction
- Alloy structure optimization

## Comparison with Other Codes
- **vs USPEX**: ParetoCSP multi-objective, ML-accelerated
- **vs CrySPY**: Different GA approach
- **vs CALYPSO**: ParetoCSP Pareto-based
- **Unique strength**: Age-fitness Pareto GA, M3GNet integration

## Community and Support
- Open-source (GitHub)
- Academic development (USC)
- Published methodology

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/sadmanomee/ParetoCSP
2. Publication: J. Mater. Inform. 4, 2 (2024)
3. arXiv: https://arxiv.org/abs/2309.06710

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE (GitHub)
- Source: OPEN
- Development: ACTIVE
- Applications: Crystal structure prediction, ML-accelerated CSP
