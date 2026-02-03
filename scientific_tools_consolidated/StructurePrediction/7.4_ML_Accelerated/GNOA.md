# GNOA

## Overview
GNOA (Graph Network + Optimization Algorithm) is a machine learning approach for crystal structure prediction that combines graph networks for energy prediction with optimization algorithms for structure search.

## Theoretical Basis
- Graph neural networks (GNN)
- Formation enthalpy prediction
- Bayesian optimization
- Particle swarm optimization
- Structure-energy correlation

## Key Capabilities
- ML-accelerated structure prediction
- Graph network energy model
- Multiple optimization algorithms
- Fast energy evaluation
- Efficient structure search

**Sources**: Nature Communications 13, 1492 (2022)

## Key Strengths

### Graph Networks:
- Structure-energy correlation
- Fast evaluation
- Transferable model

### Optimization:
- Bayesian optimization
- Particle swarm
- Efficient search

### Performance:
- Nature Communications
- Well-validated
- Competitive results

## Inputs & Outputs
- **Input formats**: Chemical composition
- **Output data types**: Predicted structures, energies

## Interfaces & Ecosystem
- **ML**: Graph neural networks
- **Optimization**: BO, PSO
- **Validation**: DFT codes

## Workflow and Usage
1. Train GNN on formation enthalpies
2. Define composition
3. Run optimization (BO/PSO)
4. Generate candidate structures
5. Validate with DFT

## Performance Characteristics
- Fast GNN evaluation
- Efficient optimization
- Good for binary systems

## Computational Cost
- GNN: fast
- Optimization: moderate
- Validation: DFT

## Best Practices
- Use appropriate GNN model
- Choose optimization algorithm
- Validate predictions
- Check structural validity

## Limitations & Known Constraints
- GNN training data dependent
- Complex systems challenging
- Requires validation

## Application Areas
- Crystal structure prediction
- Materials discovery
- High-throughput screening
- Energy prediction

## Comparison with Other Codes
- **vs CrySPY**: Different ML approach
- **vs ParetoCSP**: Different optimization
- **Unique strength**: GNN + optimization combination, Nature Comms

## Community and Support
- Academic development (Wan-Jian Yin group)
- Published methodology
- Research code

## Verification & Sources
**Primary sources**:
1. Publication: Nature Communications 13, 1492 (2022)
2. Group: http://www.comates.group/

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Documentation: AVAILABLE (paper)
- Development: ACADEMIC
- Applications: ML-accelerated CSP
