# CrySPR

## Overview
CrySPR (Crystal Structure Pre-Relaxation and PRediction) is a Python interface for implementing crystal structure pre-relaxation and prediction using machine-learning interatomic potentials (ML-IAPs).

## Theoretical Basis
- Structure generation (PyXtal)
- ML interatomic potentials (M3GNet, CHGNet, MACE)
- Structure pre-relaxation
- Random search / PSO
- Energy ranking

## Key Capabilities
- ML-IAP pre-relaxation
- PyXtal structure generation
- Multiple ML potentials supported
- Random search and PSO
- Fast screening

**Sources**: ChemRxiv (2024)

## Key Strengths

### ML Integration:
- M3GNet, CHGNet, MACE
- Fast relaxation
- DFT-level accuracy

### Workflow:
- Automated pipeline
- PyXtal integration
- Multiple search methods

### Flexibility:
- Multiple ML-IAPs
- Customizable workflow
- Python-based

## Inputs & Outputs
- **Input formats**: Chemical composition, space group (optional)
- **Output data types**: Relaxed structures, energies

## Interfaces & Ecosystem
- **Generation**: PyXtal
- **ML Potentials**: M3GNet, CHGNet, MACE
- **ASE**: Calculator interface

## Workflow and Usage
1. Define composition
2. Generate structures (PyXtal)
3. Pre-relax with ML-IAP
4. Run global search (RS/PSO)
5. Rank and validate

## Performance Characteristics
- Fast ML relaxation
- Efficient screening
- Good pre-relaxation

## Computational Cost
- ML relaxation: fast
- Search: moderate
- DFT validation: expensive

## Best Practices
- Choose appropriate ML-IAP
- Use pre-relaxation
- Validate top structures
- Compare multiple potentials

## Limitations & Known Constraints
- ML-IAP accuracy limits
- Training data dependent
- Requires validation

## Application Areas
- Crystal structure prediction
- Pre-relaxation screening
- Materials discovery
- ML-accelerated CSP

## Comparison with Other Codes
- **vs PyMCSP**: Similar purpose
- **vs CrySPY**: CrySPR pre-relaxation focused
- **Unique strength**: ML-IAP pre-relaxation interface

## Community and Support
- Open-source (GitHub)
- Recent development
- ChemRxiv preprint

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/Tosykie/CrySPR
2. ChemRxiv: https://chemrxiv.org/engage/chemrxiv/article-details/66b308a501103d79c5fd9b91

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE (GitHub)
- Source: OPEN
- Development: RECENT
- Applications: ML-accelerated CSP
