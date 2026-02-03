# MGAC / MGAC2

## Overview
MGAC (Modified Genetic Algorithm for Crystals) is a genetic algorithm package for molecular crystal structure prediction. MGAC2 is the updated version with improved algorithms.

## Theoretical Basis
- Genetic algorithm optimization
- Molecular crystal handling
- Conformational flexibility
- Space group exploration
- Force field evaluation

## Key Capabilities
- Molecular crystal CSP
- Any space group
- Multiple molecules per cell
- Conformational flexibility
- Parallel execution

**Sources**: J. Chem. Phys. 116, 5984 (2002), Chem. Phys. Lett. 623, 69 (2015)

## Key Strengths

### Molecular Crystals:
- Designed for organic systems
- Conformational flexibility
- Multiple Z' values

### Flexibility:
- Any space group
- Various force fields
- Customizable operators

### Track Record:
- Long development history
- Many publications
- Validated methodology

## Inputs & Outputs
- **Input formats**: Molecular geometry, GA parameters
- **Output data types**: Crystal structures, energies, population history

## Interfaces & Ecosystem
- **Force fields**: Various empirical potentials
- **DFT**: Optional refinement
- **Analysis**: Structure comparison

## Workflow and Usage
1. Prepare molecular geometry
2. Configure GA parameters
3. Select force field
4. Run genetic algorithm
5. Refine top structures

## Performance Characteristics
- Force field-limited
- Good parallel scaling
- Efficient for molecular crystals

## Computational Cost
- Force field: fast
- DFT refinement: expensive
- Parallelizable

## Best Practices
- Use appropriate force field
- Enable conformational search
- Refine with DFT
- Compare with experiment

## Limitations & Known Constraints
- Force field accuracy
- Molecular crystals only
- Academic distribution

## Application Areas
- Pharmaceutical polymorph prediction
- Organic crystal engineering
- Molecular materials design
- CSP blind tests

## Comparison with Other Codes
- **vs GAtor**: Similar purpose, different implementation
- **vs Genarris**: MGAC includes optimization
- **Unique strength**: Long history, molecular crystal focus

## Community and Support
- Academic development (Utah)
- Published methodology
- Available upon request

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/MGAC-group/MGAC2
2. Publications: J. Chem. Phys. 116, 5984 (2002)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE (GitHub)
- Source: AVAILABLE
- Development: MAINTAINED
- Applications: Molecular crystal CSP
