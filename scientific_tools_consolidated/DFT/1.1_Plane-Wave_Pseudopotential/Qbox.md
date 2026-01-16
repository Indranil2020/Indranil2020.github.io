# Qbox

## Official Resources
- Homepage: http://qboxcode.org/
- Documentation: http://qboxcode.org/doc/html/
- Source Repository: https://github.com/qboxcode/qbox-public
- License: GNU General Public License v2.0

## Overview
Qbox is a scalable parallel implementation of first-principles molecular dynamics based on the plane-wave, pseudopotential formalism. Developed by François Gygi at UC Davis, Qbox is specifically designed for exceptional parallel scalability on high-performance computing systems, with demonstrated efficiency on tens of thousands of processors. It excels at large-scale ab initio molecular dynamics simulations of complex systems, particularly for studying materials under extreme conditions.

**Scientific domain**: First-principles MD, plane-wave DFT, extreme scalability, HPC  
**Target user community**: Materials scientists, HPC researchers, extreme conditions, large-scale MD

## Theoretical Methods
- Kohn-Sham DFT (LDA, GGA, hybrid functionals)
- Plane-wave basis with pseudopotentials
- Norm-conserving pseudopotentials
- Born-Oppenheimer molecular dynamics
- Car-Parrinello-like dynamics
- Variable-cell dynamics
- Path integral molecular dynamics (PIMD)
- van der Waals corrections (vdW-DF)
- DFT+U for correlated systems
- Meta-GGA functionals
- Exact exchange (HFX) for hybrids
- Spin-orbit coupling
- Non-collinear magnetism
- Ultrasoft pseudopotentials (limited)

## Capabilities (CRITICAL)
- Ground-state electronic structure calculations
- Born-Oppenheimer molecular dynamics
- First-principles molecular dynamics at scale
- Geometry optimization
- Real-time time-dependent DFT
- Optical absorption spectra via RT-TDDFT
- Electronic stopping power calculations
- Forces and stress tensors
- Parallel execution on 1000s of processors
- XML-based checkpointing and I/O
- Wave function extrapolation for MD
- Variable cell dynamics
- Thermostats and barostats
- Optimized for large systems (1000+ atoms)

**Sources**: Official Qbox documentation, cited in 6/7 source lists

## Inputs & Outputs
- **Input formats**:
  - XML input files (native format)
  - Interactive command-line interface
  - Coordinate files (XYZ)
  - Pseudopotential files (XML format)
  
- **Output data types**:
  - XML output with all computed properties
  - Sample files (snapshots during MD)
  - Wavefunction checkpoints
  - Trajectory files
  - Energy and force outputs

## Interfaces & Ecosystem
- **Framework integrations**:
  - Can be interfaced with workflow systems
  - XML I/O allows custom parsing
  
- **HPC optimization**:
  - Massively parallel (MPI)
  - Optimized for Cray, IBM, Intel architectures
  - ScaLAPACK for linear algebra
  - FFTW for fast Fourier transforms
  
- **Analysis**:
  - qbox_xyz - trajectory extraction
  - Custom XML parsers
  - Standard molecular dynamics analysis tools

## Limitations & Known Constraints
- **Specialized focus**: Optimized for MD; fewer features than general DFT codes
- **Input format**: XML-based, requires learning curve
- **Pseudopotentials**: Limited to specific formats (norm-conserving primarily)
- **Post-processing**: Fewer established analysis tools compared to VASP/QE
- **Documentation**: Good but less extensive than major codes
- **Community**: Smaller user base
- **Hybrid functionals**: Limited implementation
- **Features**: Focuses on MD; fewer property calculations than general codes
- **Basis sets**: Plane-wave only
- **Platform support**: Primarily for HPC systems; not optimized for desktop use

## Computational Cost
- **Scalability**: Scaling to 100,000+ cores demonstrated (Blue Gene/Q).
- **Efficiency**: Optimized for large-scale AIMD (1000+ atoms).
- **Overhead**: Minimum system size recommended is ~50 atoms due to parallel overhead structure.

## Comparison with Other Codes
- **vs Quantum ESPRESSO**: Qbox is designed specifically for "First Principles MD at scale". QE is a general purpose swiss-army knife. Qbox is simpler but faster for its specific niche (large AIMD).
- **vs CPMD**: Qbox uses Born-Oppenheimer MD (mostly), CPMD uses Car-Parrinello. Qbox scales better on modern massive supercomputers.

## Best Practices
- **Input**: Learn the XML-based input or use the interactive command line for steering simulations.
- **HPC**: Use the `row_m`, `col_m` matrix distributions to map the calculation to the specific torus/interconnect topology of your supercomputer.
- **Restart**: Qbox's XML restart files are robust; use them frequently to checkpoint long MD runs.

## Community and Support
- **Mailing List**: `qbox-users` list available.
- **Documentation**: Online HTML manual is the primary resource.

## Verification & Sources
**Primary sources**:
1. Official website: http://qboxcode.org/
2. Documentation: http://qboxcode.org/doc/
3. GitHub repository: https://github.com/qboxcode/qbox-public
4. F. Gygi, IBM J. Res. Dev. 52, 137 (2008) - Qbox architecture
5. E. Draeger et al., J. Parallel Distrib. Comput. 106, 205 (2017) - Scaling study

**Secondary sources**:
1. Qbox tutorials and examples
2. Published AIMD applications using Qbox
3. HPC benchmarking studies
4. Confirmed in 6/7 source lists (claude, g, gr, k, m, q)

**Confidence**: CONFIRMED - Appears in 6 of 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: ACCESSIBLE
- Source code: OPEN (GitHub)
- Community support: Active (mailing list, development team)
- Academic citations: >200 (main papers)
- HPC pedigree: Developed at LLNL, optimized for supercomputers
