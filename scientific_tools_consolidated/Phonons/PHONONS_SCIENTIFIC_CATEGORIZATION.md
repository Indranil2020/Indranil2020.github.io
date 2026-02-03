# Scientific Categorization of Phonon Computational Tools

**Total Tools**: 42 verified phonon codes  
**Categorization Basis**: Theoretical methodology, physics domain, and computational approach  
**Date**: 2025-01-17

---

## 1. HARMONIC LATTICE DYNAMICS (Core Phonon Calculations)

### 1.1 Finite Displacement Method
**Theoretical Basis**: Direct force constant extraction from atomic displacements
- **Phonopy** - Industry standard, supercell finite displacement
- **PHON** - Parlinski's original finite displacement code
- **PHONON** - Wolf's implementation
- **FROPHO** - Frozen phonon by Togo
- **ASE-phonons** - ASE framework implementation
- **API_Phonons** - Python API wrapper
- **Phonopy-API** - Programmatic Phonopy access

### 1.2 Density Functional Perturbation Theory (DFPT)
**Theoretical Basis**: Linear response theory, no supercell needed
- **EPW** - Quantum ESPRESSO electron-phonon
- **PERTURBO** - Beyond EPW, advanced e-ph coupling
- **Phoebe** - Modern DFPT-based transport

### 1.3 Machine Learning Potentials for Phonons
**Theoretical Basis**: ML-accelerated force constant generation
- **hiPhive** - Force constant extraction from MD
- **mace_phonopy** - MACE ML potential bridge
- **CRYSTALpytools** - CRYSTAL code ML infrastructure

---

## 2. ANHARMONIC PHONON PHYSICS

### 2.1 Third-Order Force Constants (3-Phonon Scattering)
**Theoretical Basis**: Cubic anharmonicity, three-phonon processes
- **phono3py** - Supercell 3rd-order FCs
- **ShengBTE** - Boltzmann transport equation solver
- **thirdorder.py** - ShengBTE displacement generator
- **ALAMODE** - Compressive sensing for FCs
- **ALM** - ALAMODE force constant module

### 2.2 Fourth-Order and Higher Anharmonicity
**Theoretical Basis**: Quartic and higher-order phonon interactions
- **FourPhonon** - 4-phonon scattering processes
- **ALATDYN** - Higher-order force constants

### 2.3 DFPT-Based Anharmonic Methods
**Theoretical Basis**: Perturbation theory for anharmonicity
- **D3Q** - Quantum ESPRESSO 3rd derivatives
- **THERMAL2** - Variational BTE with D3Q

---

## 3. TEMPERATURE-DEPENDENT PHONONS

### 3.1 Temperature-Dependent Effective Potentials (TDEP)
**Theoretical Basis**: Finite-temperature renormalized force constants from AIMD
- **TDEP** - Temperature-dependent effective potential
- **a-TDEP** - Abinit TDEP implementation
- **autoplex** - Automated MLIP workflow with TDEP

### 3.2 Self-Consistent Phonon Theory
**Theoretical Basis**: Variational free energy minimization
- **SSCHA** - Stochastic self-consistent harmonic approximation

### 3.3 Molecular Dynamics-Based Methods
**Theoretical Basis**: Spectral analysis of MD trajectories
- **DynaPhoPy** - Phonon quasiparticle from MD
- **phonon-sed** - Spectral energy density from LAMMPS

---

## 4. THERMAL TRANSPORT CALCULATIONS

### 4.1 Boltzmann Transport Equation (BTE) Solvers
**Theoretical Basis**: Phonon BTE with scattering rates
- **ShengBTE** - Relaxation time approximation
- **almaBTE** - Monte Carlo BTE solver
- **phono3py** - Iterative BTE solution
- **THERMAL2** - Variational BTE
- **kALDo** - Anharmonic lattice dynamics
- **GPU_PBTE** - GPU-accelerated BTE
- **Pheasy** - Educational BTE solver
- **THERMACOND** - Thermal conductivity calculator
- **OpenBTE** - Open-source BTE

### 4.2 Molecular Dynamics Thermal Transport
**Theoretical Basis**: Green-Kubo or direct MD methods
- **GPUMD** - GPU molecular dynamics with NEP
- **Simphony** - Phonon transport simulator

### 4.3 Monte Carlo and Particle Methods
**Theoretical Basis**: Phonon particle transport
- **freepaths** - Monte Carlo phonon trajectories

---

## 5. ELECTRON-PHONON COUPLING

### 5.1 Superconductivity and Transport
**Theoretical Basis**: Electron-phonon matrix elements, Eliashberg theory
- **EPW** - Electron-phonon Wannier interpolation
- **PERTURBO** - Ab initio e-ph transport
- **epiq** - Electron-phonon interaction quantities
- **elphbolt** - Coupled electron-phonon BTE
- **ElectronPhononCoupling** - Abinit EPC analysis

### 5.2 Spectroscopy and Optical Properties
**Theoretical Basis**: Phonon-assisted optical processes
- **Phonopy-Spectroscopy** - IR and Raman spectra
- **IR2PW** - Irreducible representation analysis

### 5.3 Ultrafast Dynamics
**Theoretical Basis**: Non-equilibrium electron-phonon dynamics
- **XTANT-3** - X-ray induced thermal/non-thermal transitions
- **USER-EPH** - LAMMPS two-temperature model

---

## 6. QUASI-HARMONIC APPROXIMATION (QHA)

**Theoretical Basis**: Volume-dependent phonons, thermal expansion
- **qha** - Quasi-harmonic approximation tool
- **ATAT** - Alloy thermodynamics with phonons
- **Phonopy** - Built-in QHA module

---

## 7. ADVANCED ANALYSIS AND POST-PROCESSING

### 7.1 Symmetry and Group Theory
**Theoretical Basis**: Irreducible representations, selection rules
- **PhononIrep** - Phonon irrep analysis
- **IR2PW** - Character table analysis

### 7.2 Visualization Tools
**Theoretical Basis**: Graphical representation of phonon modes
- **phononwebsite** - Web-based 3D visualization
- **Phonopy_VESTA** - VESTA format export
- **Phonon-Vibration-Viewer** - Primitive atom visualization
- **phononplotter** - Band structure plotting
- **AtomicContributions** - Mode character analysis

### 7.3 Band Unfolding
**Theoretical Basis**: Supercell to primitive cell projection
- **easyunfold** - Spectral weight unfolding
- **KPROJ** - Electronic and phononic unfolding

### 7.4 Extended Analysis
**Theoretical Basis**: Advanced post-processing workflows
- **Phono3py-Power-Tools** - Extended Phono3py analysis
- **pwtools** - QE phonon post-processing
- **elphmod** - Electron-phonon model analysis
- **latticeDynamics** - Rigid ion model tools

---

## 8. GREEN'S FUNCTION AND NEGF METHODS

**Theoretical Basis**: Non-equilibrium Green's function formalism
- **AGF-phonon-transport** - Atomistic Green's function (MATLAB)
- **NEGF-phonon-1D-matlab** - Educational 1D NEGF

---

## 9. SPECIALIZED APPLICATIONS

### 9.1 Scattering and Disorder
**Theoretical Basis**: Phonon scattering in disordered systems
- **SCAILD** - Self-consistent anharmonic lattice dynamics
- **QSCAILD** - Quantum SCAILD

### 9.2 Educational and Prototyping
**Theoretical Basis**: Simplified implementations for learning
- **PhonTS** - Phonon transport simulator
- **NEGF-phonon-1D-matlab** - 1D NEGF tutorial
- **Pheasy** - Easy-to-use BTE solver

---

## SCIENTIFIC METHODOLOGY SUMMARY

### By Computational Approach:
1. **Supercell Methods** (18 tools): Phonopy, phono3py, ShengBTE, etc.
2. **DFPT/Linear Response** (5 tools): EPW, PERTURBO, D3Q, etc.
3. **Molecular Dynamics** (6 tools): TDEP, DynaPhoPy, GPUMD, etc.
4. **Machine Learning** (4 tools): hiPhive, autoplex, mace_phonopy, GPUMD
5. **Green's Function** (2 tools): AGF, NEGF-phonon-1D

### By Physics Domain:
1. **Harmonic Phonons** (17 tools)
2. **Anharmonic Transport** (15 tools)
3. **Electron-Phonon** (7 tools)
4. **Temperature-Dependent** (5 tools)
5. **Visualization/Analysis** (10 tools)

### By Theoretical Complexity:
1. **Harmonic Approximation**: Phonopy, PHON, ASE-phonons
2. **Quasi-Harmonic**: qha, ATAT
3. **Cubic Anharmonicity**: phono3py, ShengBTE, ALAMODE
4. **Quartic+**: FourPhonon, ALATDYN
5. **Full Anharmonic**: TDEP, SSCHA, DynaPhoPy

---

## INTEGRATION ECOSYSTEM

### DFT Code Backends:
- **VASP**: Phonopy, phono3py, ShengBTE, EPW
- **Quantum ESPRESSO**: EPW, D3Q, THERMAL2, pwtools
- **Abinit**: PERTURBO, a-TDEP, ElectronPhononCoupling
- **CASTEP**: Phonopy, Phono3py-Power-Tools
- **CRYSTAL**: CRYSTALpytools

### Framework Integration:
- **ASE**: ASE-phonons, multiple interfaces
- **pymatgen**: AtomicContributions, easyunfold
- **LAMMPS**: phonon-sed, USER-EPH, GPUMD
- **Wannier90**: EPW, PERTURBO

---

## VERIFICATION STATUS

- **Fully Verified**: 40/42 tools
- **Uncertain**: 2/42 (DMDW, RTDW - likely typos or internal codes)
- **Active Development**: 35/42 tools
- **Educational/Tutorial**: 5/42 tools

---

**Notes**:
- Many tools overlap in capabilities (e.g., Phonopy can do QHA, phono3py can do harmonic)
- Categorization based on primary/unique scientific contribution
- Some tools are modules/extensions of larger packages
- ML-based methods are emerging as a major trend (2020+)
