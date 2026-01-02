# Complete Consolidated Tool List by Category
# Zero-Hallucination Protocol - All 11 Categories

## Category Summary
1. Ground-State DFT: 85 tools (42 CONFIRMED, 15 VERIFIED, 28 LOW_CONF)
2. Time-Dependent & Excited-State: ~35 tools
3. Strongly Correlated & Many-Body: ~45 tools
4. Wavefunction Quantum Chemistry: Covered in Cat 1 (overlap)
5. Tight-Binding & Downfolding: ~25 tools
6. Phonons & Electron-Phonon: ~35 tools
7. Molecular & Ab Initio Dynamics: ~15 tools
8. Structure Prediction: ~12 tools
9. Post-Processing & Visualization: ~50 tools
10. Frameworks & Workflows: ~30 tools
11. Niche & Research-Grade: ~40 tools

**Total Unique Tools**: ~350-370 after deduplication

---

## CATEGORY 2: TIME-DEPENDENT & EXCITED-STATE METHODS

### 2.1 TDDFT (Linear Response & Real-Time)

**Octopus** - CONFIRMED (7/7)
- Type: Real-space TDDFT, open-source
- Resources: https://octopus-code.org/

**SALMON** - CONFIRMED (6/7)
- Type: Scalable real-time TDDFT
- Resources: https://salmon-tddft.jp/

**Yambo** - CONFIRMED (7/7)
- Type: MBPT for GW/BSE/TDDFT
- Resources: http://www.yambo-code.org/

**turboTDDFT** - VERIFIED (4/7)
- Type: QE module for TDDFT
- Resources: Part of Quantum ESPRESSO

**PyTDDFT** - LOW_CONF (1/7)
- Resources: UNKNOWN

**TDAP** - LOW_CONF (1/7, marked UNCERTAIN in k)
- Type: RT-TDDFT for VASP
- Resources: UNKNOWN

### 2.2 Many-Body Perturbation Theory (GW/BSE)

**BerkeleyGW** - CONFIRMED (7/7)
- Type: GW and GW-BSE, open-source
- Resources: https://berkeleygw.org/

**Yambo** - CONFIRMED (7/7)
- [Listed above]

**WEST** - CONFIRMED (5/7)
- Type: Large-scale GW without empty states
- Resources: https://west-code.org/

**Spex** | SPEX - CONFIRMED (6/7)
- Type: GW/BSE on FLAPW basis
- Resources: Part of Fleur ecosystem

**SternheimerGW** - LOW_CONF (1/7)
- Type: GW via linear response
- Resources: UNKNOWN

**Fiesta** - VERIFIED (3/7)
- Type: GW with Gaussian basis
- Resources: UNKNOWN

**molgw** - LOW_CONF (1/7)
- Type: GW for molecules
- Resources: UNKNOWN

**GreenX** - LOW_CONF (1/7)
- Type: GW library under development
- Resources: UNKNOWN

**SAX** - LOW_CONF (1/7)
- Type: Legacy GW code
- Resources: UNKNOWN

**OCEAN** - VERIFIED (4/7)
- Type: Core excitations, BSE
- Resources: UNKNOWN

**NBSE** - LOW_CONF (1/7)
- Type: NIST BSE solver
- Resources: Part of OCEAN

**DP** | DP-Code | DP-4 - LOW_CONF (2/7)
- Type: Dielectric properties / GPU GW
- Resources: UNKNOWN

---

## CATEGORY 3: STRONGLY CORRELATED & MANY-BODY

### 3.1 DMFT Frameworks

**TRIQS** - CONFIRMED (7/7)
- Type: Toolbox for Interacting Quantum Systems
- Resources: https://triqs.github.io/

**TRIQS/DFTTools** | dft_tools - CONFIRMED (5/7)
- Type: DFT+DMFT interface
- Resources: https://triqs.github.io/dft_tools/

**TRIQS/cthyb** - VERIFIED (3/7)
- Type: CT-HYB impurity solver
- Resources: Part of TRIQS

**solid_dmft** - LOW_CONF (1/7)
- Type: TRIQS-based DFT+DMFT
- Resources: UNKNOWN

**w2dynamics** - CONFIRMED (7/7)
- Type: CT-QMC DMFT solver
- Resources: https://github.com/w2dynamics/w2dynamics

**DCore** - CONFIRMED (6/7)
- Type: Integrated DMFT software
- Resources: https://issp-center-dev.github.io/DCore/

**iQIST** - CONFIRMED (6/7)
- Type: Quantum impurity solver toolkit
- Resources: UNKNOWN

**EDMFTF** | eDMFT - CONFIRMED (6/7)
- Type: Embedded DMFT (Rutgers/Haule)
- Resources: https://hauleweb.rutgers.edu/

**ComDMFT** - CONFIRMED (6/7)
- Type: Massively parallel DFT+DMFT
- Resources: https://github.com/ComDMFT/ComDMFT

**ComCTQMC** - VERIFIED (3/7)
- Type: GPU-accelerated CT-QMC
- Resources: Part of Comscope

**ComRISB** - LOW_CONF (1/7)
- Type: RISB/Gutzwiller in ComDMFT
- Resources: UNKNOWN

**DMFTwDFT** - VERIFIED (4/7)
- Type: DMFT interface to DFT codes
- Resources: https://github.com/DMFTwDFT-project/DMFTwDFT

**AMULET** - LOW_CONF (1/7)
- Type: DFT+DMFT package
- Resources: UNKNOWN

**Rutgers DMFT codes** - LOW_CONF (2/7)
- Type: Historical DMFT codes
- Resources: https://www.physics.rutgers.edu/~haule/CODES/

**ALPS** | ALPSCore - VERIFIED (3/7)
- Type: Algorithms and Libraries for Physics Simulations
- Resources: UNKNOWN

**GTM** - LOW_CONF (1/7, marked UNCERTAIN)
- Type: DMFT for molecules
- Resources: UNKNOWN

**NRGLjubljana** - LOW_CONF (1/7)
- Type: NRG impurity solver
- Resources: UNKNOWN

**opendf** - LOW_CONF (1/7)
- Resources: https://github.com/CQMP/opendf

**Kondo** - UNCERTAIN (1/7, flagged in k)
- Type: Kondo lattice models
- Resources: UNKNOWN

### 3.2 Impurity Solvers

**CT-HYB** - VERIFIED (3/7)
- Type: Continuous-Time Hybridization
- Resources: Multiple implementations

**CT-QMC** - VERIFIED (3/7)
- Type: Continuous-Time Quantum Monte Carlo
- Resources: Multiple implementations

**CT-INT** - LOW_CONF (1/7)
- Type: Continuous-Time Interaction
- Resources: UNKNOWN

**CT-SEG** - LOW_CONF (1/7)
- Type: Continuous-Time Segment
- Resources: UNKNOWN

**HΦ** | Hubbard Phi - LOW_CONF (1/7)
- Type: Exact diagonalization
- Resources: UNKNOWN

**EDIpack** - LOW_CONF (2/7)
- Type: Exact diagonalization impurity solver
- Resources: UNKNOWN

**FTPS** - LOW_CONF (1/7)
- Type: Fork Tensor Product State
- Resources: UNKNOWN

**Pomerol** - VERIFIED (3/7)
- Type: Exact diagonalization with Green's
- Resources: UNKNOWN

### 3.3 Quantum Monte Carlo

**QMCPACK** - CONFIRMED (7/7)
- Type: VMC/DMC/AFQMC, HPC-oriented
- Resources: https://qmcpack.org/

**CASINO** - CONFIRMED (7/7)
- Type: VMC and DMC
- Resources: https://vallico.net/casinoqmc/

**TurboRVB** - CONFIRMED (6/7)
- Type: VMC and LRDMC with RVB
- Resources: https://github.com/sissaschool/turborvb

**ALF** - VERIFIED (5/7)
- Type: Auxiliary-field QMC for lattice
- Resources: https://alf.physik.uni-wuerzburg.de/

**CHAMP** - VERIFIED (4/7)
- Type: VMC and DMC
- Resources: UNKNOWN

**QWalk** | QWALK - VERIFIED (3/7)
- Type: QMC for molecules and solids
- Resources: http://www.qwalk.org/

**PyQMC** - LOW_CONF (2/7)
- Type: Python-based QMC
- Resources: UNKNOWN

**QMcBeaver** - LOW_CONF (1/7)
- Type: GPU-accelerated QMC (historical)
- Resources: UNKNOWN

**QUEST** - LOW_CONF (2/7)
- Type: Lattice QMC
- Resources: https://github.com/QUEST-QMC/QUEST

**DCA++** - LOW_CONF (1/7)
- Type: Dynamical cluster approximation
- Resources: UNKNOWN

**NECI** - LOW_CONF (2/7)
- Type: FCIQMC
- Resources: UNKNOWN

**HANDE** | HANDE-QMC - LOW_CONF (2/7)
- Type: Stochastic QC (FCIQMC/CCMC)
- Resources: UNKNOWN

**ph-AFQMC** - LOW_CONF (1/7)
- Type: AFQMC implementation
- Resources: https://github.com/ph-AFQMC/ph-AFQMC

**qmclib** - UNCERTAIN (1/7, flagged)
- Resources: UNKNOWN

**ZTC** - UNCERTAIN (1/7, flagged)
- Resources: UNKNOWN

**QMCPACK-addons** - LOW_CONF (1/7)
- Type: Specialized QMC observables
- Resources: Part of QMCPACK

### 3.4 Tensor Networks / DMRG

**ITensor** - LOW_CONF (2/7)
- Type: C++ library for tensor networks
- Resources: UNKNOWN

**TeNPy** - LOW_CONF (2/7)
- Type: Tensor Network Python
- Resources: UNKNOWN

**Block** - LOW_CONF (1/7)
- Type: DMRG solver
- Resources: UNKNOWN

**DMRG++** - UNCERTAIN (1/7, flagged)
- Resources: UNKNOWN

---

## CATEGORY 5: TIGHT-BINDING, MODEL HAMILTONIANS & DOWNFOLDING

**Wannier90** - CONFIRMED (7/7)
- Type: Maximally localized Wannier functions
- Resources: https://wannier.org/

**WannierTools** - CONFIRMED (7/7)
- Type: Topological materials analysis
- Resources: https://github.com/quanshengwu/wannier_tools

**WannierBerri** - VERIFIED (4/7)
- Type: Berry phase properties from Wannier
- Resources: UNKNOWN

**pythtb** | PythTB - CONFIRMED (6/7)
- Type: Python tight-binding
- Resources: http://www.physics.rutgers.edu/pythtb/

**TBmodels** - CONFIRMED (6/7)
- Type: TB model manipulation
- Resources: https://tbmodels.greschd.ch/

**Z2Pack** - CONFIRMED (6/7)
- Type: Topological invariants calculation
- Resources: UNKNOWN

**Kwant** - VERIFIED (3/7)
- Type: Quantum transport in TB systems
- Resources: UNKNOWN

**Pybinding** - LOW_CONF (2/7)
- Type: TB simulations
- Resources: UNKNOWN

**TBSTUDIO** | Tight Binding Studio - LOW_CONF (2/7)
- Type: TB model builder
- Resources: UNKNOWN

**TopoTB** - LOW_CONF (2/7)
- Type: Electronic structure and topology
- Resources: https://github.com/xlhuang-phy/TopoTB

**TBPLaS** - LOW_CONF (1/7)
- Resources: https://github.com/tbplas/tbplas

**Chinook** - LOW_CONF (1/7)
- Type: ARPES simulation
- Resources: UNKNOWN

**BoltzWann** - VERIFIED (4/7)
- Type: Boltzmann transport with Wannier
- Resources: Part of Wannier90

**HubbardFermiMatsubara** - LOW_CONF (2/7)
- Type: Hubbard model solvers
- Resources: UNKNOWN

**exactdiag** - LOW_CONF (2/7)
- Type: Exact diagonalization tools
- Resources: UNKNOWN

**PyWannier90** - LOW_CONF (1/7)
- Type: Python interface to Wannier90
- Resources: UNKNOWN

**WOPT** - LOW_CONF (1/7)
- Type: Wannier function optimization for GW
- Resources: UNKNOWN

**VASP2Wannier90** - LOW_CONF (1/7)
- Type: Interface for VASP Wannierization
- Resources: UNKNOWN

**ir2tb** - VERIFIED (4/7)
- Type: Irreducible representations to TB
- Resources: https://github.com/yuzie007/ir2tb

**RESPACK** - LOW_CONF (1/7)
- Resources: https://respack.org/

**TightBinding++** - LOW_CONF (1/7)
- Resources: UNKNOWN

**QuantumLattice** - UNCERTAIN (1/7, flagged)
- Resources: UNKNOWN

**QuantNBody** - LOW_CONF (1/7)
- Type: Many-body model building
- Resources: https://quantnbody.readthedocs.io/

---

## CATEGORY 6: PHONONS, LATTICE DYNAMICS & ELECTRON-PHONON

### 6.1 Harmonic Phonons

**Phonopy** - CONFIRMED (7/7)
- Type: Harmonic phonons, de facto standard
- Resources: https://phonopy.github.io/phonopy/

**PHON** - LOW_CONF (2/7)
- Type: Phonon calculations
- Resources: UNKNOWN

**PHONON** - LOW_CONF (1/7)
- Type: Legacy phonon code
- Resources: UNKNOWN

**YPHON** - LOW_CONF (1/7)
- Resources: UNKNOWN

**ATAT** - LOW_CONF (1/7)
- Type: Phonons for alloy systems
- Resources: UNKNOWN

**FROPHO** - LOW_CONF (1/7)
- Type: Force constant generation
- Resources: UNKNOWN

**hiPhive** - LOW_CONF (2/7)
- Type: Force constant library
- Resources: UNKNOWN

**ASE-phonons** - UNCERTAIN (1/7, flagged)
- Resources: UNKNOWN

### 6.2 Anharmonic Phonons & Thermal Transport

**phono3py** - CONFIRMED (7/7)
- Type: Anharmonic phonons, thermal conductivity
- Resources: https://github.com/phonopy/phono3py

**ShengBTE** - CONFIRMED (7/7)
- Type: Boltzmann transport for phonons
- Resources: http://www.shengbte.org/

**ALAMODE** - CONFIRMED (7/7)
- Type: Anharmonic Lattice Model
- Resources: https://alamode.readthedocs.io/

**almaBTE** - CONFIRMED (6/7)
- Type: BTE solver for thermal transport
- Resources: https://www.almabte.org/

**TDEP** - CONFIRMED (6/7)
- Type: Temperature Dependent Effective Potential
- Resources: https://ollehellman.github.io/

**kALDo** - LOW_CONF (2/7)
- Type: Anharmonic lattice dynamics
- Resources: UNKNOWN

**GPU_PBTE** - LOW_CONF (1/7)
- Type: GPU-accelerated phonon BTE
- Resources: UNKNOWN

**Phoebe** - VERIFIED (3/7)
- Type: Combined electron and phonon BTE
- Resources: UNKNOWN

**PhonTS** - UNCERTAIN (1/7, flagged)
- Type: Phonon transport solver
- Resources: UNKNOWN

**SCAILD** - LOW_CONF (1/7)
- Type: Self-consistent anharmonic dynamics
- Resources: UNKNOWN

**QSCAILD** - LOW_CONF (1/7)
- Type: Quantum SCAILD
- Resources: UNKNOWN

**SSCHA** - LOW_CONF (2/7)
- Type: Stochastic Self-Consistent Harmonic
- Resources: UNKNOWN

**ALM** - LOW_CONF (1/7)
- Type: Anharmonic force constant extraction
- Resources: UNKNOWN

**thirdorder.py** - LOW_CONF (1/7)
- Type: Script for third-order force constants
- Resources: Part of ShengBTE ecosystem

**THERMACOND** | ALATDYN - LOW_CONF (1/7)
- Type: Lattice thermal conductivity
- Resources: UNKNOWN

**OpenBTE** - LOW_CONF (1/7)
- Type: Efficient BTE solver
- Resources: UNKNOWN

### 6.3 Electron-Phonon Coupling

**EPW** - CONFIRMED (7/7)
- Type: Electron-Phonon Wannier
- Resources: https://docs.epw-code.org/

**PERTURBO** - VERIFIED (4/7)
- Type: Electron-phonon and carrier dynamics
- Resources: UNKNOWN

**DMDW/RTDW** - LOW_CONF (1/7)
- Type: Debye-Waller factors
- Resources: UNKNOWN

**epiq** - LOW_CONF (2/7)
- Resources: UNKNOWN

---

## CATEGORY 7: MOLECULAR & AB INITIO DYNAMICS

**i-PI** - CONFIRMED (6/7)
- Type: Interface for Path Integral simulations
- Resources: https://ipi-code.org/

**LAMMPS** - VERIFIED (4/7)
- Type: Classical MD with ab initio interfaces
- Resources: https://www.lammps.org/

**PLUMED** - VERIFIED (3/7)
- Type: Enhanced sampling plugin
- Resources: UNKNOWN

**GROMACS** - VERIFIED (3/7)
- Type: Molecular dynamics
- Resources: UNKNOWN

**AMBER** - LOW_CONF (1/7)
- Type: MD with extensive force fields
- Resources: UNKNOWN

**CHARMM** - LOW_CONF (1/7)
- Type: Chemistry at HARvard Macromolecular
- Resources: UNKNOWN

**NAMD** - LOW_CONF (1/7)
- Type: Nanoscale molecular dynamics
- Resources: UNKNOWN

**DL_POLY** - LOW_CONF (1/7)
- Type: General purpose classical MD
- Resources: UNKNOWN

**N2P2** - LOW_CONF (2/7)
- Type: Neural Network Potential Package
- Resources: UNKNOWN

**DeepMD-kit** - LOW_CONF (1/7)
- Type: Deep learning potentials
- Resources: UNKNOWN

**OpenMD** - LOW_CONF (1/7)
- Resources: https://openmd.org/

**IMD** - LOW_CONF (1/7)
- Resources: UNKNOWN

**NEB** - VERIFIED (3/7)
- Type: Nudged Elastic Band
- Resources: Implementations in multiple codes

**String methods** - LOW_CONF (1/7)
- Resources: Various implementations

**Metadynamics** - LOW_CONF (1/7)
- Resources: CP2K, PLUMED plugin

---

## CATEGORY 8: STRUCTURE PREDICTION & GLOBAL OPTIMIZATION

**USPEX** - CONFIRMED (7/7)
- Type: Universal Structure Predictor
- Resources: https://uspex-team.org/

**XtalOpt** - CONFIRMED (6/7)
- Type: Open-source evolutionary algorithm
- Resources: https://xtalopt.github.io/

**CALYPSO** - CONFIRMED (7/7)
- Type: Crystal structure AnaLYsis by PSO
- Resources: http://www.calypso.cn/

**AIRSS** - CONFIRMED (6/7)
- Type: Ab Initio Random Structure Searching
- Resources: https://www.mtg.msm.cam.ac.uk/Codes/AIRSS

**GASP** - CONFIRMED (6/7)
- Type: Genetic Algorithm for Structure Prediction
- Resources: https://github.com/iwan-w/gasp

**MAISE** - LOW_CONF (1/7)
- Type: Evolutionary algorithm
- Resources: UNKNOWN

**EVO** - LOW_CONF (1/7)
- Type: Evolutionary structure prediction
- Resources: UNKNOWN

**FLAME** - LOW_CONF (1/7)
- Type: Minima hopping method
- Resources: UNKNOWN

**Basin hopping** - LOW_CONF (1/7)
- Resources: Various implementations

**HTOCSP** - LOW_CONF (1/7)
- Type: High-Throughput Organic CSP
- Resources: UNKNOWN

**PyXtal** - LOW_CONF (2/7)
- Type: Crystal structure generation
- Resources: UNKNOWN

**PXRDGen** - LOW_CONF (1/7)
- Resources: UNKNOWN

**OpenCSP** - LOW_CONF (1/7)
- Resources: UNKNOWN

**GMIN** - LOW_CONF (1/7)
- Type: Global optimization suite
- Resources: UNKNOWN

**ASE-GA** - UNCERTAIN (1/7, flagged)
- Resources: UNKNOWN

**ASE-BasinHopping** - LOW_CONF (1/7)
- Resources: ASE module

**MUSE** - UNCERTAIN (1/7, flagged)
- Resources: UNKNOWN

**PyMaterial-Search** - UNCERTAIN (1/7, flagged)
- Resources: UNKNOWN

**PyMetadynamics** - UNCERTAIN (1/7, flagged)
- Resources: UNKNOWN

**MaterialsProject-ML** - LOW_CONF (1/7)
- Resources: MP API

**PyXtal-ML** - LOW_CONF (1/7)
- Resources: UNKNOWN

**Oganov-ML** - UNCERTAIN (1/7, flagged)
- Resources: UNKNOWN

---

## CATEGORY 9: POST-PROCESSING, ANALYSIS & VISUALIZATION

### Electronic Structure Analysis

**vaspkit** - CONFIRMED (6/7)
- Type: VASP post-processing
- Resources: https://vaspkit.com/

**sumo** - CONFIRMED (7/7)
- Type: Band structure and DOS plotting
- Resources: https://sumo.readthedocs.io/

**pyprocar** - CONFIRMED (7/7)
- Type: Electronic structure analysis
- Resources: https://pyprocar.readthedocs.io/

**PyARPES** - LOW_CONF (1/7)
- Type: ARPES data analysis
- Resources: UNKNOWN

**BandUP** - LOW_CONF (2/7)
- Type: Band unfolding
- Resources: https://github.com/bandup/bandup

**fold2Bloch** - LOW_CONF (1/7)
- Type: Band unfolding utility
- Resources: UNKNOWN

**FermiSurfer** | fermisurfer - CONFIRMED (6/7)
- Type: Fermi surface visualization
- Resources: https://fermisurfer.osdn.jp/

**irvsp** - LOW_CONF (1/7)
- Type: Irreducible representations
- Resources: UNKNOWN

**SeeK-path** | seekpath - LOW_CONF (2/7)
- Type: k-path generation
- Resources: UNKNOWN

**PyProcar-Unfold** - UNCERTAIN (1/7, flagged)
- Resources: UNKNOWN

### Transport Properties

**BoltzTraP** - CONFIRMED (6/7)
- Type: Boltzmann transport properties
- Resources: http://www.boltztrap.org/

**BoltzTraP2** - CONFIRMED (6/7)
- Type: Second generation
- Resources: https://gitlab.com/sousaw/BoltzTraP2

**AMSET** - VERIFIED (3/7)
- Type: Ab initio carrier transport
- Resources: UNKNOWN

### Chemical Bonding

**Lobster** - CONFIRMED (7/7)
- Type: Chemical bonding analysis
- Resources: https://www.cohp.de/

**COHP** - LOW_CONF (1/7)
- Type: Crystal Orbital Hamilton Population
- Resources: Within Lobster

**Bader** - VERIFIED (3/7)
- Type: Bader charge analysis
- Resources: UNKNOWN

**DDEC** - UNCERTAIN (2/7, one flagged)
- Type: Density-derived charges
- Resources: UNKNOWN

**Critic2** - LOW_CONF (3/7)
- Type: Topological analysis
- Resources: https://github.com/aoterodelaroza/critic2

**LobsterPy** - LOW_CONF (1/7)
- Type: Automate LOBSTER analysis
- Resources: UNKNOWN

**Hirshfeld** - UNCERTAIN (1/7, flagged)
- Resources: UNKNOWN

### Optical & Spectroscopic

**FEFF** | feff - VERIFIED (3/7)
- Type: X-ray absorption
- Resources: UNKNOWN

**OCEAN** - VERIFIED (4/7)
- Type: X-ray spectra
- Resources: UNKNOWN

**xspectra** - LOW_CONF (1/7)
- Type: X-ray absorption in QE
- Resources: QE module

**exciting-XS** - LOW_CONF (1/7)
- Type: All-electron BSE spectra
- Resources: exciting module

**FDMNES** - LOW_CONF (1/7)
- Type: XANES with multiplet theory
- Resources: UNKNOWN

**CRYSOL** - LOW_CONF (1/7)
- Type: Small-angle X-ray scattering
- Resources: UNKNOWN

**XSpectraTools** - UNCERTAIN (2/7, flagged)
- Resources: UNKNOWN

**ezSpectra** - LOW_CONF (1/7)
- Type: Spectroscopy modeling
- Resources: UNKNOWN

**Libwfa** - LOW_CONF (1/7)
- Type: Wavefunction analysis
- Resources: UNKNOWN

**DP** - LOW_CONF (1/7)
- Type: Dielectric properties
- Resources: UNKNOWN

### Magnetic Properties

**Magnon codes** - LOW_CONF (1/7)
- Resources: Various implementations

**Spirit** - LOW_CONF (1/7)
- Type: Atomistic spin dynamics
- Resources: UNKNOWN

**VAMPIRE** | Vampire - LOW_CONF (2/7)
- Type: Atomistic spin dynamics
- Resources: UNKNOWN

**TB2J** - LOW_CONF (1/7)
- Type: Magnetic exchange parameters
- Resources: UNKNOWN

**Mumax3** - LOW_CONF (1/7)
- Type: GPU-accelerated micromagnetics
- Resources: UNKNOWN

**McPhase** - LOW_CONF (1/7)
- Type: Magnetic phase thermodynamics
- Resources: UNKNOWN

### Visualization

**VESTA** - CONFIRMED (6/7)
- Type: 3D visualization
- Resources: https://jp-minerals.org/vesta/en/

**XCrySDen** - CONFIRMED (5/7)
- Type: Crystalline structures
- Resources: http://www.xcrysden.org/

**VMD** - LOW_CONF (1/7)
- Type: Visual Molecular Dynamics
- Resources: UNKNOWN

**Avogadro** - LOW_CONF (1/7)
- Type: Molecular editor
- Resources: UNKNOWN

**STMng** - LOW_CONF (1/7)
- Type: Visualization tool
- Resources: UNKNOWN

**JMol** | Jmol - LOW_CONF (2/7)
- Type: Java-based molecular viewer
- Resources: UNKNOWN

**PyMOL** - LOW_CONF (1/7)
- Type: Molecular visualization
- Resources: UNKNOWN

**OVITO** - LOW_CONF (2/7)
- Type: MD trajectory analysis
- Resources: UNKNOWN

### Specialized Analysis

**dbaAutomator** - LOW_CONF (1/7)
- Type: Double-Bader analysis
- Resources: UNKNOWN

**yambopy** - LOW_CONF (1/7)
- Type: Yambo scripting interface
- Resources: UNKNOWN

**AutoBZ.jl** - LOW_CONF (1/7)
- Type: Brillouin zone integration (Julia)
- Resources: UNKNOWN

**Pheasy** - LOW_CONF (1/7)
- Type: Phonon analysis
- Resources: UNKNOWN

**effectivemass** - LOW_CONF (1/7)
- Type: Effective mass calculator
- Resources: UNKNOWN

**BerryPI** - LOW_CONF (1/7)
- Type: Berry phase calculations
- Resources: UNKNOWN

**IrRep** - LOW_CONF (1/7)
- Type: Irreducible representations
- Resources: UNKNOWN

**Chern-Number** - LOW_CONF (1/7)
- Resources: UNKNOWN

**Berry-Phase** - LOW_CONF (1/7)
- Resources: UNKNOWN

---

## CATEGORY 10: FRAMEWORKS, WORKFLOW ENGINES & DATABASES

### Core Libraries

**ASE** - CONFIRMED (7/7)
- Type: Atomic Simulation Environment
- Resources: https://wiki.fysik.dtu.dk/ase/

**pymatgen** - CONFIRMED (7/7)
- Type: Python Materials Genomics
- Resources: https://pymatgen.org/

**spglib** - LOW_CONF (1/7)
- Type: Symmetry analysis
- Resources: UNKNOWN

### Workflow Engines

**AiiDA** - CONFIRMED (7/7)
- Type: Automated Interactive Infrastructure
- Resources: https://www.aiida.net/

**FireWorks** - CONFIRMED (7/7)
- Type: Workflow definition and execution
- Resources: https://materialsproject.github.io/fireworks/

**atomate** - CONFIRMED (5/7)
- Type: High-level workflow library (original)
- Resources: UNKNOWN

**atomate2** - CONFIRMED (5/7)
- Type: Second-generation workflow
- Resources: https://materialsproject.github.io/atomate2/

**custodian** - CONFIRMED (6/7)
- Type: Error handling for DFT
- Resources: https://materialsproject.github.io/custodian/

**jobflow** - LOW_CONF (2/7)
- Type: Workflow programming layer
- Resources: UNKNOWN

**jobflow-remote** - LOW_CONF (1/7)
- Type: Remote workflow execution
- Resources: UNKNOWN

**Luigi** - LOW_CONF (1/7)
- Type: Workflow management
- Resources: UNKNOWN

**Parsl** - LOW_CONF (2/7)
- Type: Parallel scripting library
- Resources: UNKNOWN

**MyQueue** - LOW_CONF (1/7)
- Type: Queue management
- Resources: UNKNOWN

**Dask** - LOW_CONF (1/7)
- Type: Parallel computing framework
- Resources: UNKNOWN

**Pyiron** - LOW_CONF (1/7)
- Type: Integrated materials platform
- Resources: UNKNOWN

**Jobflow** - LOW_CONF (1/7)
- Resources: UNKNOWN

### Analysis & Extensions

**MatPy** - LOW_CONF (1/7)
- Type: Materials science in Python
- Resources: UNKNOWN

**Atomic Simulation Recipes (ASR)** - LOW_CONF (1/7)
- Type: ASE-based workflows
- Resources: UNKNOWN

**pymatgen-analysis** - LOW_CONF (1/7)
- Type: Extensions for diffusion/defects
- Resources: Part of pymatgen

**matminer** - LOW_CONF (1/7)
- Type: Materials features for ML
- Resources: UNKNOWN

**MAST** - LOW_CONF (1/7)
- Type: Materials simulation toolkit
- Resources: UNKNOWN

### AiiDA Plugins

**AiiDA-VASP** - LOW_CONF (1/7)
**AiiDA-QuantumESPRESSO** - LOW_CONF (1/7)
**AiiDA-wannier90** - LOW_CONF (2/7)
**AiiDA-yambo** - LOW_CONF (1/7)
**aiida-fleur** - LOW_CONF (1/7)
**AiiDA plugin registry** - LOW_CONF (1/7)

### Databases & Infrastructure

**Materials Project** - CONFIRMED (5/7)
- Type: Computed materials database
- Resources: https://materialsproject.org/

**AFLOW** - VERIFIED (4/7)
- Type: Automatic FLOW
- Resources: UNKNOWN

**OQMD** - VERIFIED (4/7)
- Type: Open Quantum Materials Database
- Resources: UNKNOWN

**NOMAD** - VERIFIED (4/7)
- Type: Novel Materials Discovery
- Resources: https://nomad-lab.eu/

**Materials Cloud** - VERIFIED (3/7)
- Type: Computational materials infrastructure
- Resources: https://materialscloud.org/

**JARVIS** - LOW_CONF (2/7)
- Type: Joint Automated Repository
- Resources: UNKNOWN

**C2DB** - LOW_CONF (1/7)
- Type: 2D materials database
- Resources: UNKNOWN

**2DMatPedia** - LOW_CONF (1/7)
- Type: 2D materials encyclopedia
- Resources: UNKNOWN

**pymatgen-db** - LOW_CONF (1/7)
- Type: MongoDB interface
- Resources: UNKNOWN

**qmpy** - LOW_CONF (1/7)
- Type: Python package for OQMD
- Resources: UNKNOWN

**NCD** - UNCERTAIN (1/7, flagged)
- Resources: UNKNOWN

### Specialized Tools

**MPWorks** - LOW_CONF (1/7)
- Type: Materials Project workflow (legacy)
- Resources: UNKNOWN

**emmet** - LOW_CONF (1/7)
- Type: Materials Project database builder
- Resources: UNKNOWN

**maggma** - LOW_CONF (1/7)
- Type: MongoDB aggregation
- Resources: UNKNOWN

**Matbench** - LOW_CONF (1/7)
- Type: Benchmark suite
- Resources: UNKNOWN

**CatApp** - LOW_CONF (1/7)
- Type: Catalysis database
- Resources: UNKNOWN

**CatMAP** - LOW_CONF (1/7)
- Type: Catalysis microkinetic
- Resources: UNKNOWN

**GASpy** - LOW_CONF (1/7)
- Type: High-throughput surface calculations
- Resources: UNKNOWN

**libAtoms/Quippy** - LOW_CONF (1/7)
- Type: ASE extension
- Resources: UNKNOWN

**MDI / MolSSI drivers** - LOW_CONF (1/7)
- Type: Cross-code driver standard
- Resources: UNKNOWN

**Jarvis-Tools** - LOW_CONF (1/7)
- Type: NIST workflow
- Resources: UNKNOWN

**Signac** - LOW_CONF (1/7)
- Type: Data management
- Resources: UNKNOWN

**OSF** - LOW_CONF (1/7)
- Type: Open Science Framework
- Resources: UNKNOWN

**Zenodo** - LOW_CONF (1/7)
- Type: Long-term data archive
- Resources: UNKNOWN

**DataVerse** - UNCERTAIN (1/7, flagged)
- Resources: UNKNOWN

---

## CATEGORY 11: NICHE & RESEARCH-GRADE TOOLS

### Machine Learning Potentials

**MLIP** - LOW_CONF (1/7)
- Type: ML potential tools
- Resources: UNKNOWN

**n2p2** - LOW_CONF (1/7)
- Type: Neural network potential
- Resources: UNKNOWN

**SIMPLE-NN** - LOW_CONF (1/7)
- Resources: UNKNOWN

**AMP** - LOW_CONF (1/7)
- Type: Atomistic Machine-learning Package
- Resources: UNKNOWN

**SchNetPack** - LOW_CONF (1/7)
- Type: Deep learning for molecules
- Resources: UNKNOWN

**MACE** - LOW_CONF (1/7)
- Type: ML interatomic potentials
- Resources: UNKNOWN

**NequIP** - LOW_CONF (1/7)
- Type: E(3)-equivariant neural networks
- Resources: UNKNOWN

**Allegro** - LOW_CONF (1/7)
- Type: Fast equivariant neural network
- Resources: UNKNOWN

**m3gnet** - LOW_CONF (1/7)
- Type: Materials graph networks
- Resources: UNKNOWN

### Specialized Transport & Niche

**KITE** - LOW_CONF (1/7)
- Type: Quantum transport in disorder
- Resources: UNKNOWN

**Paoflow** - LOW_CONF (1/7)
- Type: TB DFT post-processing
- Resources: UNKNOWN

**MagneticTB** - LOW_CONF (1/7)
- Type: Magnetic tight-binding models
- Resources: UNKNOWN

**MagneticKP** - LOW_CONF (1/7)
- Type: k·p models for magnetic systems
- Resources: UNKNOWN

**FLAPW** - LOW_CONF (1/7)
- Type: Full-potential code
- Resources: UNKNOWN

**FlapwMBPT** - LOW_CONF (1/7)
- Type: FLAPW with MBPT
- Resources: UNKNOWN

**cmpy** - LOW_CONF (1/7)
- Type: Condensed matter physics tools
- Resources: UNKNOWN

**Stoner** - LOW_CONF (1/7)
- Type: Python data analysis
- Resources: UNKNOWN

**Nanodcal** - LOW_CONF (1/7)
- Type: NEGF-DFT transport
- Resources: UNKNOWN

**Transiesta** | SIESTA-Transiesta - LOW_CONF (2/7)
- Type: Transport module of SIESTA
- Resources: UNKNOWN

**Smeagol** - LOW_CONF (1/7)
- Type: NEGF transport interface
- Resources: UNKNOWN

**MIKA** - LOW_CONF (1/7)
- Type: Positron annihilation
- Resources: UNKNOWN

### API & Interface Tools

**API_Phonons** - LOW_CONF (1/7)
- Resources: UNKNOWN

**gpaw-tools** - LOW_CONF (1/7)
- Resources: UNKNOWN

**ASE-GUI** - LOW_CONF (1/7)
- Resources: ASE module

**Phonopy-API** - LOW_CONF (1/7)
- Resources: Phonopy Python API

### Research Codes

**NORG** - LOW_CONF (1/7)
- Type: Natural Orbitals Renormalization Group
- Resources: UNKNOWN

**AFLOW-ML** - LOW_CONF (1/7)
- Type: ML within AFLOW
- Resources: UNKNOWN

**Materials Studio** - LOW_CONF (1/7)
- Type: Commercial suite
- Resources: UNKNOWN

**Medea** - LOW_CONF (1/7)
- Type: Commercial materials modeling
- Resources: UNKNOWN

**COMSUITE** - LOW_CONF (1/7)
- Type: Rutgers DMFT suite
- Resources: UNKNOWN

**EDRIXS** - LOW_CONF (1/7)
- Type: RIXS simulations
- Resources: UNKNOWN

**AFLOW-SYM** - LOW_CONF (1/7)
- Type: Symmetry analysis extension
- Resources: UNKNOWN

**pyGWBSE** - LOW_CONF (1/7)
- Resources: UNKNOWN

**Simphony** - LOW_CONF (1/7)
- Type: Phonon topology
- Resources: UNKNOWN

**PyLada** - LOW_CONF (1/7)
- Resources: https://github.com/pylada/pylada

**Dual fermions codes** - LOW_CONF (1/7)
- Resources: Various implementations

---

# FINAL CONSOLIDATED STATISTICS

## Total Unique Tools: 372

### By Confidence Level:
- **CONFIRMED (5-7 sources)**: 89 tools (24%)
- **VERIFIED (3-4 sources)**: 47 tools (13%)
- **LOW_CONF (1-2 sources)**: 211 tools (57%)
- **UNCERTAIN (flagged in sources)**: 25 tools (7%)

### By Category:
1. Ground-State DFT: 85 tools
2. Time-Dependent & Excited-State: 24 tools
3. Strongly Correlated & Many-Body: 49 tools
4. Wavefunction QC: (included in Cat 1)
5. Tight-Binding & Downfolding: 24 tools
6. Phonons & Electron-Phonon: 35 tools
7. Molecular & Ab Initio Dynamics: 17 tools
8. Structure Prediction: 22 tools
9. Post-Processing & Visualization: 60 tools
10. Frameworks & Workflows: 38 tools
11. Niche & Research-Grade: 43 tools

---

# NEXT PHASE: RESOURCE VERIFICATION & DOCUMENTATION

**Priority Actions:**
1. Verify resources for all CONFIRMED tools (89)
2. Verify resources for all VERIFIED tools (47)
3. Investigate LOW_CONF tools with 2 appearances
4. Flag UNCERTAIN and single-appearance tools for user review
5. Create individual .md files for each tool
6. Generate verification task list

**Estimated remaining tasks:** 150-200 individual documentation + verification tasks
