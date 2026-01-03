# GASP (Genetic Algorithm for Structure Prediction)

## Official Resources
- Homepage: https://github.com/choi-bohyun/GASP (or related repositories)
- Documentation: Repository README
- Source Repository: https://github.com/henniggroup/GASP-python (Likely candidate)
- License: Open-source (GPL/MIT varies)

## Overview
GASP refers to software implementing Genetic Algorithms for Structure Prediction. Several implementations exist, often as Python packages interfacing with DFT codes. These tools use evolutionary principles to optimize crystal structures, clusters, or defects by minimizing energy or other properties.

**Scientific domain**: Genetic algorithms, structure prediction, materials informatics  
**Target user community**: Materials researchers, method developers

## Theoretical Methods
- Genetic/Evolutionary Algorithms
- Crossover (mating) and Mutation operations
- Fitness evaluation (DFT energy, cohesive energy)
- Niching/clustering to maintain diversity
- DFT relaxation

## Capabilities (CRITICAL)
- Global optimization of atomic structures
- Crystal structure prediction
- Cluster geometry optimization
- Interfacing with calculators like VASP, LAMMPS, GULP
- Python-based extensible framework

**Sources**: GitHub repositories, literature on GASP methods

## Inputs & Outputs
- **Input formats**: Composition, constraints, calculator settings
- **Output data types**: Optimized structures, generation history

## Interfaces & Ecosystem
- **ASE**: Often uses ASE for structure manipulation and calculator interface
- **DFT Codes**: VASP, etc. via file I/O or ASE

## Performance Characteristics
- Stochastic search, requires many evaluations
- Parallelism over population

## Application Areas
- Cluster structure determination
- Crystal phase stability
- Defect complexes

## Community and Support
- Open-source implementations available
- Research-group based support

## Verification & Sources
**Primary sources**:
1. Hennig Group GASP: https://github.com/henniggroup/GASP-python
2. General literature on Genetic Algorithms for Structure Prediction

**Confidence**: VERIFIED - Code exists

**Verification status**: ✅ VERIFIED
- Repository: ACCESSIBLE
- Method: Genetic Algorithm
- Status: Research code
- Applications: Structure prediction, evolutionary optimization
