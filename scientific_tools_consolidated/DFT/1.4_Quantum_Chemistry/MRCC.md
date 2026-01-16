# MRCC

## Official Resources
- Homepage: https://www.mrcc.hu/
- Documentation: https://www.mrcc.hu/manual/
- Source Repository: Available to licensees
- License: Academic and commercial licenses available

## Overview
MRCC (Multi-Reference Coupled Cluster) is a specialized quantum chemistry program suite featuring arbitrary-order coupled cluster and configuration interaction methods. Developed by Mihály Kállay and collaborators at Budapest University of Technology and Economics, it is renowned for implementing the highest-order correlation methods available, including fully automated arbitrary-order CC and CI implementations generated via string-based equations.

**Scientific domain**: High-order correlation methods, multireference systems, benchmark calculations  
**Target user community**: Researchers requiring very high accuracy or exploring high-order correlation methods

## Theoretical Methods
- Hartree-Fock (RHF, UHF, ROHF)
- Density Functional Theory (DFT with Libxc)
- Møller-Plesset perturbation theory (MP2 through arbitrary order)
- Coupled Cluster up to arbitrary order:
  - CCSD, CCSD(T), CC(T), CCSDT
  - CCSDTQ, CCSDTQP, CCSDTQPH
  - CC(n), CC[n] series
- Multireference CC (Mk-MRCC, BW-MRCC)
- Configuration Interaction up to full CI
- Linear Response CC for excited states
- Equation-of-Motion CC (EOM-CCSD, EOM-CCSDT)
- Symmetry-Adapted Perturbation Theory (SAPT)
- Local Natural Orbital methods (LNO-CCSD(T))
- Explicitly correlated F12 methods
- Analytical gradients (CCSD, CCSD(T), LNO-CCSD(T))
- Relativistic methods (DKH, X2C)
- Automated generation of CC equations

## Capabilities (CRITICAL)
- Ground-state electronic structure
- Very high-order coupled cluster (up to CCSDTQPH)
- Automated coupled cluster equation generation
- String-based many-body theory
- Geometry optimization with analytic gradients
- Excited states via EOM-CC at various orders
- Local correlation for large molecules (LNO-CCSD(T))
- Explicitly correlated F12 methods
- Multi-reference CC calculations
- Molecular properties
- Benchmark-quality calculations
- Interface to external programs for integrals
- Efficient parallelization (OpenMP + MPI)
- Automated focal-point analysis
- Composite thermochemistry protocols

**Sources**: Official MRCC documentation, cited in 6/7 source lists

## Key Strengths

### Arbitrary-Order CC:
- Automated equation generation
- String-based many-body theory
- CC up to CCSDTQPH (sextuple excitations)
- CI up to full CI
- Systematic hierarchy

### Local Correlation (LNO):
- Local Natural Orbital framework
- LNO-CCSD(T) for large molecules
- Near-linear scaling
- Analytical gradients
- Production quality

### Benchmark Accuracy:
- Highest-order correlation methods
- Full CI limit approachable
- Thermochemical protocols
- Reference calculations
- Method validation

### Automated Methods:
- Focal-point analysis
- Composite methods
- Basis set extrapolation
- Automated workflows
- Error estimation

## Inputs & Outputs
- **Input formats**:
  - MRCC input file (MINP)
  - XYZ coordinate files
  - Interface inputs from CFOUR, Molpro, ORCA
  
- **Output data types**:
  - Detailed output files
  - Energies, gradients
  - Correlation energies by order
  - Property calculations
  - Wavefunction analysis

## Interfaces & Ecosystem
- **Integral interfaces**:
  - Built-in integral code (default)
  - CFOUR for advanced integrals
  - Molpro interface
  - ORCA interface (>v5.0)
  - Dirac interface for relativistic integrals
  - PSI4 interface
  
- **Standalone capabilities**:
  - Complete standalone operation
  - Built-in SCF, DFT
  - Built-in basis sets
  
- **Utilities**:
  - dmrcc - main driver
  - Automated protocol execution
  - Focal-point automation


## Workflow and Usage

### Input Format (MINP):
MRCC uses a keyword-based `MINP` input file.

```
calc=CCSD
basis=cc-pVDZ
mem=1000MB

geom=xyz
3

O 0.0 0.0 0.0
H 0.0 0.7 0.0
H 0.7 0.0 0.0
```

### Running MRCC:
```bash
dmrcc > minp.out
```

### Common Tasks:
- **Single Point**: `calc=CCSD(T)`
- **Optimization**: `geom=opt`
- **LNO Method**: `localcc=on`
- **F12**: `densfit=on` and F12 keywords

## Advanced Features

### Arbitrary-Order CC:
- `calc=CC(n)` allows easy access to high orders (e.g., `calc=CCSDT`)
- Automated derivation of equations
- Validated against FCI for small systems
- Up to sextuple excitations (or more)

### LNO-CCSD(T):
- Linear Scaling Local Natural Orbital CCSD(T)
- "Black-box" accuracy control (Tight/Normal/Loose)
- Applicable to systems with 100+ atoms
- Analytical gradients for geometry optimization

### F12 Methods:
- Compatible with arbitrary order CC
- Drastically reduces basis set error
- Approach FCI limit with manageable basis sets

### Relativistic Effects:
- `calc=X2C` or `calc=DKH2`
- Essential for heavy element accuracy
- Compatible with LNO methods

### Automated Protocols:
- `calc=HEAT` for thermochemistry
- `calc=W4`
- Automated basis set extrapolation
- Composite energy schemes

## Performance Characteristics
- **Accuracy**: The "Reference" code for high-order coupled cluster
- **Speed**: Generally slower than specialized low-order codes (PSI4/Molpro/ORCA) for standard CCSD(T), but unique for high orders
- **Scalability**: Hybrid MPI/OpenMP parallelization; LNO methods scale linearly
- **Memory**: Arbitrary order methods scale factorially in memory/disk
- **Disk I/O**: Very heavy for high-order methods

## Computational Cost
- **CCSD(T)**: O(N^7), standard
- **CCSDT**: O(N^8)
- **CCSDTQ**: O(N^10)
- **LNO-CCSD(T)**: Linear scaling, O(N), breaks even ~30-50 atoms
- **Arbitrary Order**: Grows extremely fast, limited to <10 atoms for very high orders

## Comparison with Other Codes
- **vs CFOUR**: CFOUR better for analytical derivatives/properties; MRCC better for higher Order energies and local correlation.
- **vs Molpro**: Molpro is "gold standard" for MRCI; MRCC is "gold standard" for high-order CC.
- **vs ORCA**: ORCA's DLPNO is similar to MRCC's LNO; MRCC offers higher order canonical CC benchmarks.
- **Unique strength**: Arbitrary-order Coupled Cluster (CCSDT, CCSDTQ, ...), LNO-CCSD(T) for large molecules.

## Best Practices

### High-Order Calculations:
- Use small basis sets first to test feasibility
- Estimate memory requirements carefully
- Use `restart` options for long jobs

### LNO-CCSD(T):
- Use `lc_ortho=tight` for benchmark accuracy
- Check domain sizes
- Use density fitting (`densfit=on`) for speed

### Methods:
- Use `calc=CCSD(T)` for standard chemical accuracy
- Use `calc=CCSDT` only for benchmarking small systems
- Use F12 basis sets with F12 methods

## Community and Support
- **Support**: Active email support from developers
- **Manual**: Detailed description of keywords
- **Development**: Kállay group (Budapest)
- **License**: Academic (free/low cost) and Commercial

## Application Areas

### Benchmark Calculations:
- Reference energies
- Method validation
- Convergence studies
- Correlation hierarchy
- FCI extrapolation

### High-Accuracy Thermochemistry:
- Atomization energies
- Reaction barriers
- Heats of formation
- Isomerization energies
- Sub-kJ/mol accuracy

### Large Molecule Applications:
- LNO-CCSD(T) for 50-100 atoms
- Organic molecules
- Biomolecule fragments
- Drug-like molecules
- Noncovalent interactions

### Method Development:
- High-order CC research
- Perturbation theory studies
- Local correlation development
- Multi-reference theory

## Limitations & Known Constraints
- **Registration required**: Free for academics but requires registration
- **Molecular focus**: Not designed for periodic systems
- **System size**: High-order CC limited to very small molecules (<10 atoms)
- **Memory**: Very high-level methods extremely memory-intensive
- **Computational cost**: High-order CC scales steeply (factorial-like)
- **Basis sets**: Gaussian-type; large basis required for accuracy
- **Learning curve**: Steep; requires expert knowledge
- **Documentation**: Good but assumes high-level theory knowledge
- **Parallelization**: Efficient OpenMP+MPI hybrid
- **Platform**: Linux, macOS, Windows

## Verification & Sources
**Primary sources**:
1. Official website: https://www.mrcc.hu/
2. Manual: https://www.mrcc.hu/manual/
3. M. Kállay et al., J. Chem. Phys. 152, 074107 (2020) - MRCC program
4. M. Kállay and P. R. Surján, J. Chem. Phys. 115, 2945 (2001) - Arbitrary-order CC
5. P. R. Nagy, M. Kállay, J. Chem. Phys. 150, 104101 (2019) - LNO-CCSD(T)

**Secondary sources**:
1. MRCC manual and examples
2. Published benchmark calculations
3. High-accuracy thermochemistry studies
4. Confirmed in 6/7 source lists (claude, g, gr, k, m, q)

**Confidence**: CONFIRMED - Appears in 6 of 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: ACCESSIBLE
- Software: Free for academics (registration required)
- Community support: Active (email support)
- Academic citations: >1,000
- Unique capability: Automated arbitrary-order CC equations, highest-order correlation, LNO local correlation
- Specialized strength: Arbitrary-order CC (up to CCSDTQPH), string-based many-body theory, local correlation (LNO-CCSD(T)), benchmark thermochemistry, automated protocols
