# WEST

## Official Resources
- Homepage: https://west-code.org/
- Documentation: https://west-code.org/doc/
- Source Repository: https://github.com/west-code-development/West
- License: GNU General Public License v3.0

## Overview
WEST (Without Empty STates) is a massively parallel software for large-scale electronic structure calculations within many-body perturbation theory. It implements GW approximation and TDDFT for computing quasiparticle energies and optical spectra without the need for empty states, enabling efficient calculations for large systems using density functional perturbation theory techniques.

**Scientific domain**: Many-body perturbation theory, GW calculations, optical properties, large-scale systems  
**Target user community**: Researchers studying electronic excitations, band structures, and optical properties of materials

## Theoretical Methods
- GW approximation (G₀W₀, eigenvalue self-consistent GW)
- Perturbative GW without empty states
- Density functional perturbation theory (DFPT) formalism
- Lanczos algorithm for frequency integration
- Contour deformation techniques
- TDDFT linear response (planned)
- Bethe-Salpeter equation (under development)
- Static COHSEX approximation
- Plasmon-pole models
- Full-frequency integration

## Capabilities (CRITICAL)
- GW quasiparticle energies and corrections
- Band structure corrections (gaps, alignments)
- Density of states with many-body effects
- Large-scale GW calculations (1000+ atoms)
- No empty states required (reduced computational cost)
- Frequency-dependent dielectric function
- Static and dynamic screening
- Self-energy calculations
- Spectral functions
- Interface with Quantum ESPRESSO
- Massively parallel (MPI+OpenMP)
- GPU acceleration (CUDA)
- Memory-efficient algorithms
- Linear scaling approaches
- Efficient dielectric matrix calculation

**Sources**: Official WEST documentation (https://west-code.org/), cited in 5/7 source lists

## Key Innovations

### No Empty States Paradigm:
- Traditional GW requires summing over empty (unoccupied) states
- WEST uses DFPT to avoid explicit empty states
- Computes response directly via Lanczos chains
- Dramatically reduces memory and computational cost
- Enables GW for systems with 1000+ atoms

### Computational Efficiency:
- Linear scaling with system size (for local excitations)
- Massively parallel implementation
- GPU acceleration for key kernels
- Reduced memory footprint
- Efficient frequency integration

### Scalability:
- Demonstrated on systems with >1000 atoms
- Excellent parallel scaling to 10,000+ cores
- Suitable for HPC and exascale computing
- Production runs on leadership facilities

## Inputs & Outputs
- **Input formats**:
  - west.in (WEST input file)
  - Quantum ESPRESSO DFT outputs
  - Charge density and wavefunctions from QE
  - Pseudopotential files
  - k-point grids from DFT
  
- **Output data types**:
  - Quasiparticle energies (eV)
  - Self-energy corrections
  - Band structures (GW-corrected)
  - Density of states
  - Dielectric function
  - Convergence data
  - JSON output files

## Interfaces & Ecosystem
- **Quantum ESPRESSO integration**:
  - Primary DFT backend (exclusive interface)
  - Reads QE wavefunctions and charge density
  - Uses QE pseudopotentials
  - Tight integration with PWscf
  
- **Analysis tools**:
  - Python scripts for post-processing (westpy)
  - Jupyter notebooks for tutorials
  - Plotting utilities
  - Data extraction tools
  
- **HPC optimization**:
  - MPI parallelization (FFT, bands, k-points)
  - OpenMP threading
  - GPU offloading (CUDA)
  - Optimized for modern architectures

## Workflow and Usage

### Typical GW Workflow:

1. **DFT Ground State**:
   - Run Quantum ESPRESSO DFT calculation
   - Converge charge density and wavefunctions
   - Generate necessary output files

2. **WEST Setup**:
   - Prepare west.in input file
   - Specify k-points for GW corrections
   - Set convergence parameters
   - Define frequency grid

3. **GW Calculation**:
   - Run wstat.x (screening calculation)
   - Run wfreq.x (self-energy and QP energies)
   - Monitor convergence

4. **Post-Processing**:
   - Extract quasiparticle energies
   - Plot corrected band structures
   - Analyze density of states
   - Compare with experiments

### Example Input:
```
input_west:
  qe_prefix: 'silicon'
  west_prefix: 'silicon_gw'
  outdir: './'

wstat_control:
  wstat_calculation: 'S'
  n_pdep_eigen: 50
  n_pdep_times: 4

wfreq_control:
  wfreq_calculation: 'XWGQ'
  n_pdep_eigen_to_use: 50
  qp_bandrange: [1,8]
```

## Advanced Features

### Lanczos Algorithm:
- Efficiently computes dielectric response
- Avoids explicit sum over empty states
- Iterative approach with controlled convergence
- Balances accuracy and computational cost

### Frequency Integration:
- Contour deformation for efficiency
- Full-frequency or plasmon-pole options
- Adaptive frequency grids
- Analytic continuation schemes

### Parallelization Strategy:
- Band/k-point parallelization
- FFT-based parallelization
- Image parallelization for independent calculations
- Hybrid MPI+OpenMP
- GPU acceleration for linear algebra

### Memory Optimization:
- Out-of-core algorithms for large systems
- Distributed memory model
- Efficient storage of response functions
- Checkpoint and restart capabilities

## Performance Characteristics
- **Typical GW calculation**: Hours to days
- **System size**: Up to 1000+ atoms demonstrated
- **Parallel efficiency**: >80% to 10,000+ cores
- **GPU speedup**: 2-5x for key operations
- **Memory**: Significantly lower than traditional GW
- **Scaling**: Near-linear for local excitations

## Computational Comparison
- **vs Traditional GW**: 10-100x faster for large systems
- **vs BerkeleyGW**: WEST better for very large systems
- **vs Yambo**: Similar capabilities, different algorithms
- **Unique strength**: No empty states, extreme scalability

## Limitations & Known Constraints
- **Quantum ESPRESSO dependency**: Requires QE for DFT part
- **Pseudopotentials**: Limited to norm-conserving
- **GW approximation**: Does not include vertex corrections
- **Learning curve**: Moderate; requires GW knowledge
- **Input complexity**: Multiple input files and convergence parameters
- **Documentation**: Good but assumes familiarity with GW theory
- **BSE**: Under development, not yet production-ready
- **k-point sampling**: Dense meshes required for accuracy
- **Frequency grids**: Convergence testing necessary
- **Platform**: Linux/Unix; HPC recommended

## Application Areas

### Materials Science:
- Band gap predictions and corrections
- Heterojunctions and interfaces
- Defect levels in semiconductors
- 2D materials and nanostructures

### Large-Scale Systems:
- Nanoparticles and clusters
- Supercells with defects
- Disordered and amorphous materials
- Biological systems (in development)

### Method Development:
- GW algorithm benchmarking
- Scaling studies
- Accuracy assessments
- Frequency integration schemes

### Comparison with Experiments:
- Photoemission spectroscopy (PES/ARPES)
- Inverse photoemission (IPES)
- Optical absorption (when BSE available)
- Electronic band structures

## Best Practices

### Convergence Testing:
- Number of Lanczos iterations (n_pdep_eigen)
- Frequency grid (n_pdep_times)
- k-point sampling
- Cutoff energies from DFT
- Check all parameters systematically

### Computational Strategy:
- Start with small test systems
- Use plasmon-pole for initial estimates
- Full-frequency for final production
- Monitor memory and time carefully
- Use checkpointing for long runs

### Validation:
- Compare with experiments
- Cross-check with other GW codes
- Test convergence thoroughly
- Verify self-consistency

## Community and Development
- Active development by Argonne National Laboratory
- Regular workshops and tutorials
- Open-source on GitHub
- Growing user community
- Integration with DOE computing facilities

## Verification & Sources
**Primary sources**:
1. Official website: https://west-code.org/
2. Documentation: https://west-code.org/doc/
3. GitHub repository: https://github.com/west-code-development/West
4. M. Govoni and G. Galli, J. Chem. Theory Comput. 11, 2680 (2015) - WEST method
5. V. Yu et al., J. Chem. Theory Comput. 16, 5035 (2020) - Large-scale GW
6. M. Govoni et al., npj Comput. Mater. 5, 96 (2019) - GPU acceleration

**Secondary sources**:
1. WEST tutorials and workshops
2. Published large-scale GW applications
3. Argonne Leadership Computing Facility resources
4. Confirmed in 5/7 source lists (claude, g, gr, k, m)

**Confidence**: CONFIRMED - Appears in 5 of 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE and ACCESSIBLE
- Source code: OPEN (GitHub, GPL v3)
- Community support: Active (mailing list, workshops, GitHub)
- Academic citations: >100
- Active development: Regular releases, DOE support
- Benchmark validation: Extensive comparisons with experiments and other GW codes
- HPC deployment: Production use on leadership facilities
