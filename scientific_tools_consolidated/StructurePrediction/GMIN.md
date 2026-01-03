# GMIN

## Official Resources
- Homepage: http://www-wales.ch.cam.ac.uk/GMIN/
- Documentation: http://www-wales.ch.cam.ac.uk/GMIN/doc/
- Source Repository: Distributed via website (Academic license)
- License: Academic License

## Overview
GMIN is a program for finding global minima of potential energy surfaces, developed by the Wales Group at the University of Cambridge. It implements the basin-hopping algorithm (and variants) to locate the lowest energy structures of clusters and molecules. It is a reference implementation for basin-hopping methods.

**Scientific domain**: Global optimization, cluster physics, energy landscapes  
**Target user community**: Physical chemists, cluster researchers

## Theoretical Methods
- Basin-Hopping (Monte Carlo Minimization)
- Parallel Tempering
- Genetic Algorithms
- Euclidean transformation coordinates
- Global optimization

## Capabilities (CRITICAL)
- Finding global minima of clusters (Lennard-Jones, Morse, etc.)
- Optimization of biomolecules (with appropriate potentials)
- Interface with various empirical potentials
- Rigorous global optimization workflows

**Sources**: Wales Group website

## Inputs & Outputs
- **Input formats**: Control file (keywords), coordinates
- **Output data types**: Lowest minima, database of structures

## Interfaces & Ecosystem
- **Potentials**: Built-in empirical potentials, interface to AMBER/CHARMM
- **OPTIM/PATHSAMPLE**: Companion codes for transition states and pathways

## Performance Characteristics
- Highly efficient for cluster optimization
- Parallelized basin-hopping

## Application Areas
- Cluster structure prediction
- Protein folding landscapes
- Molecular conformation

## Community and Support
- David Wales Group (Cambridge)
- Academic user base

## Verification & Sources
**Primary sources**:
1. Homepage: http://www-wales.ch.cam.ac.uk/GMIN/

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: AVAILABLE
- Source: ACADEMIC
- Development: ACTIVE (Wales Group)
- Applications: Global optimization, basin hopping, clusters
