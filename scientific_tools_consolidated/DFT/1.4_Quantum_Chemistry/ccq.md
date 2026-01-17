# ccq

## Official Resources
- Homepage: https://github.com/jjgoings/ccq
- Documentation: In repository
- Source Repository: https://github.com/jjgoings/ccq
- License: MIT License

## Overview
ccq is a coupled-cluster code designed for quantum chemistry calculations. It implements standard coupled-cluster methods including CCD, CCSD, CCSDT, and CCSDTQ, providing a clean implementation for research and educational purposes.

**Scientific domain**: Coupled-cluster quantum chemistry  
**Target user community**: Researchers studying coupled-cluster methods and correlation effects

## Theoretical Methods
- Coupled Cluster Doubles (CCD)
- Coupled Cluster Singles and Doubles (CCSD)
- Coupled Cluster Singles, Doubles, Triples (CCSDT)
- Coupled Cluster Singles, Doubles, Triples, Quadruples (CCSDTQ)
- Perturbative corrections
- Reference: RHF

## Capabilities (CRITICAL)
- Full CCSD implementation
- Full CCSDT implementation
- Full CCSDTQ implementation
- Clean Python code
- Tensor contraction
- Energy calculations
- Amplitude equations
- Iterative solvers
- Benchmark calculations

## Key Strengths

### Method Variety:
- Multiple CC truncations
- Full implementation (not approximate)
- High-level correlation
- Benchmark quality

### Implementation:
- Clear code structure
- Educational value
- Extensible design
- Python-based

### Research Tool:
- Method development
- Testing algorithms
- Comparison studies
- Teaching purposes

## Inputs & Outputs
- **Input formats**:
  - Molecular integrals
  - Python API
  
- **Output data types**:
  - Correlation energies
  - Amplitudes
  - Convergence data

## Interfaces & Ecosystem
- **Integral sources**: External integral packages
- **NumPy**: Array computations
- **SciPy**: Numerical methods

## Performance Characteristics
- **Speed**: Standard CC scaling
- **Accuracy**: Full CC methods
- **System size**: Small molecules
- **Memory**: CC amplitudes storage

## Computational Cost
- **CCSD**: O(N^6)
- **CCSDT**: O(N^8)
- **CCSDTQ**: O(N^10)
- **Typical**: Benchmarks on small systems

## Limitations & Known Constraints
- **System size**: Limited by CC scaling
- **Production**: Research/educational focus
- **Methods**: CC ground state
- **Large systems**: Not suitable

## Comparison with Other Codes
- **vs CFOUR/MRCC**: ccq educational, production codes optimized
- **vs ccpy**: Different implementations
- **vs PySCF CC**: ccq standalone
- **Unique strength**: Clean full CC implementations

## Application Areas

### Benchmarking:
- Reference calculations
- Method comparison
- Basis set studies

### Education:
- Learning CC theory
- Algorithm understanding
- Code modification

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/jjgoings/ccq
2. CC theory literature
3. Research applications

**Confidence**: VERIFIED
- Source code: OPEN (GitHub, MIT)
- Implementation: Standard CC
