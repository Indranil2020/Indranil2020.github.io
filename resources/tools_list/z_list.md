1. Ground-state electronic structure (DFT & variants)
Subtask: Plane‑wave / pseudopotential DFT

Quantum ESPRESSO – integrated suite of plane‑wave/pseudopotential DFT codes (pw.x, ph.x, etc.). Community: solid‑state physics and materials science. Ecosystem: standalone, with ASE/AiiDA/atomate/pymatgen interfaces.
quantum-espresso
VASP – plane‑wave PAW DFT code (commercial). Community: broad CMP/MS. Ecosystem: VASP Wiki, VASP ML‑FF into LAMMPS, pymatgen/ASE/AiiDA/atomate workflows.
vasp
+1
ABINIT – plane‑wave pseudopotential DFT code with GW/BSE and DMFT. Community: CMP/MS. Ecosystem: open source, interfaces to Phonopy, AiiDA, ASE.
abinit
+1
CASTEP – plane‑wave pseudopotential DFT (commercial/academic). Community: materials science (especially UK). Ecosystem: used in various high‑throughput and phonon workflows.
CPMD – plane‑wave/pseudopotential DFT with Car–Parrinello MD. Community: condensed matter, liquids, interfaces. Ecosystem: historically academic, now open on GitHub with MIT license.
cpmd
+1
CP2K – Gaussian & Plane Wave (GPW/GAPW) DFT and dynamics. Community: molecules, condensed matter, interfaces. Ecosystem: open source, interfaces to ASE, i‑PI, AiiDA.
cp2k
+2
PWscf – plane‑wave DFT module inside Quantum ESPRESSO (pw.x). Same communities and ecosystem as QE.
quantum-espresso
Qbox – real‑space plane‑wave/pseudopotential DFT code. Community: US CMP/MS. Ecosystem: standalone, with ASE support (listed on ASE calculators).
PARSEC – real‑space pseudopotential DFT code. Community: CMP. Ecosystem: interfaces via ASE and other I/O layers (see ASE calculators).
JDFTx – all‑electron/real‑space DFT code with focus on electrochemistry and interfaces. Community: CMP/electrochemistry. Ecosystem: open source, pymatgen I/O support.
yambo-code
Subtask: Localized basis / NAO / LCAO DFT

SIESTA – numerical atomic‑orbital DFT with order‑N scaling. Community: large‑scale solid‑state and nanostructures. Ecosystem: open source, interfaces to Phonopy, ASE, AiiDA.
FHI-aims – all‑electron NAO DFT code. Community: CMP/MS with emphasis on accuracy and large systems. Ecosystem: academic/commercial licensing; tight integration with ASE, pymatgen, AiiDA.
fhi-aims
+1
OpenMX – order‑N DFT with pseudo‑atomic orbitals. Community: large‑scale materials in Japan and globally. Ecosystem: open source (GPL), ASE OpenMX calculator, AiiDA plugins.
nsc.liu
+2
ONETEP – linear‑scaling DFT using NGWFs (localized Wannier‑like functions). Community: large‑scale CMP/MS in UK/Europe. Ecosystem: open source, supports Wannier‑like downfolding.
docs.onetep
+1
CONQUEST – linear‑scaling DFT with localized orbitals. Community: large‑scale CMP/MS. Ecosystem: open source, designed for millions of atoms, interfaced to Wannier90 workflows.
order-n
+1
BigDFT – wavelet‑based DFT. Community: CMP/MS; used within European projects. Ecosystem: open source, linked to ABINIT and other codes via wavelet libraries.
abinit
GPAW – real‑space/PAW DFT with plane‑wave and LCAO modes, part of ASE ecosystem. Community: CMP and quantum chemistry. Ecosystem: tightly integrated with ASE; used in TDDFT workflows.
wiki.fysik.dtu
DFTB+ – DFT‑based tight binding (approximate DFT) with many extensions. Community: large‑scale materials and molecules. Ecosystem: standalone or library; interfaces to ASE, i‑PI, pymatgen/ASE workflows.
dftbplus
+2
Subtask: All‑electron FP‑(L)APW / LMTO / KKR

WIEN2k – commercial FP‑LAPW code. Community: CMP, especially magnetism and heavy elements. Ecosystem: used with BoltzTraP, DMFTwWien2k, Wannier90, etc.
physics.rutgers
FLEUR – open‑source FP‑LAPW code. Community: CMP and correlated materials. Ecosystem: MaX EU ecosystem; Spex MBPT package built on FLAPW basis.
flapw
+2
Elk – open‑source FP‑LAPW code. Community: solid‑state physics. Ecosystem: open source, used with Phonopy, ASE interfaces.
elk.sourceforge
+1
exciting – FP‑(L)APW+lo code with strong focus on excitations (TDDFT, GW, e‑ph). Community: CMP, especially optical properties. Ecosystem: open source, supports Wannier90, DFPT phonons, and has TDDFT/GW modules.
exciting-code
+1
Questaal – DFT code based on LMTO basis, interfaced to DMFT via Haule’s CTQMC. Community: correlated materials. Ecosystem: open source, combined with Wannier90 for DFT+DMFT; interfaces to Wannier workflows.
SPR‑KKR – KKR Green’s function code for DFT and beyond (spin‑polarized relativistic). Community: magnetism, transport. Ecosystem: primarily standalone, used in specialized communities; see KKR review pages.
github
+1
JuKKR – Jülich KKR code (open source). Community: CMP, magnetism and disorder. Ecosystem: open source, part of JuDFT family.
github
Subtask: Real‑space FEM / high‑order FEM DFT

NESSIE – all‑electron DFT with real‑space finite‑element discretization; also real‑time TDDFT. Community: HPC and methodology. Ecosystem: designed for large‑scale parallelism; research‑grade.
nessie-code
+1
DFT‑FE – massively parallel adaptive finite‑element DFT code. Community: HPC materials science. Ecosystem: standalone HPC code for DFT ground state.
sciencedirect
Subtask: Semi‑empirical and “tight‑binding‑like” quantum methods

xtb – extended tight‑binding (GFN‑xTB) program. Community: quantum chemistry and large‑scale molecular simulations. Ecosystem: open source; integrates with ORCA and other QC workflows; available via conda.
github
+1
MOPAC – semi‑empirical quantum chemistry. Community: chemistry, large systems. Ecosystem: standalone, long‑standing code; listed in Wikipedia baseline.
AMS/DFTB (SCM ADF suite) – DFTB module within Amsterdam Modeling Suite. Community: chemistry and materials. Ecosystem: part of commercial AMS suite; used alongside ADF and BAND.
scm
Notes, cross‑check, uncertainties:

Pass 1 (Wikipedia): All plane‑wave/NAO/all‑electron codes listed above appear in the Wikipedia “List of quantum chemistry and solid‑state software” table (ABINIT, VASP, CASTEP, CPMD, CP2K, SIESTA, FHI-aims, OpenMX, ONETEP, CONQUEST, BigDFT, GPAW, WIEN2k, FLEUR, OpenMolcas/MOLCAS, etc.).
Pass 2 (major per subfield): Additional real‑space FEM (NESSIE, DFT‑FE) and KKR codes are outside or only partially covered by Wikipedia, so they are added here from primary sources.
nessie-code
+1
Pass 3 (ecosystems): ASE calculators page confirms many of these as interfaced calculators (ABINIT, CP2K, Elk, FHI-aims, OpenMX, etc.).
vasp
+1
Pass 4 (phonon/transport/topology): Phonopy, phono3py, EPW, BoltzWann, Z2Pack, WannierTools (below) explicitly rely on DFT packages above; no conflicts.
Pass 5 (small/niche): NESSIE and DFT‑FE are relatively niche HPC codes; completeness in this subclass cannot be guaranteed because such methodologically specialized codes are numerous and not always in public registries.
2. Time‑dependent & excited‑state methods
2.1 TDDFT & real‑time propagation
Octopus – real‑space grid TDDFT code with real‑time propagation, optics, strong‑field dynamics. Community: CMP and molecular TDDFT. Ecosystem: open source, supports Ehrenfest and electron dynamics.
GPAW TDDFT – real‑time TDDFT module within GPAW. Community: users of GPAW for excited states. Ecosystem: integrated in GPAW/ASE workflows.
wiki.fysik.dtu
NWChem RT‑TDDFT – real‑time TDDFT module in NWChem. Community: quantum chemistry and materials. Ecosystem: NWChem with ASE calculators; RT‑TDDFT documented in NWChem docs.
nwchemgit.github
SALMON – real‑time real‑space TDDFT code for electron dynamics in molecules and solids. Community: ultrafast spectroscopy, light–matter interaction. Ecosystem: open source (ARTED + GCEED unified).
salmon-tddft
+1
CP2K RT‑TDDFT – recent real‑time TDDFT implementation in CP2K with k‑point sampling. Community: CMP/MS and excited‑state dynamics. Ecosystem: part of CP2K, compatible with CP2K/GPW workflows.
pubs.acs
Yambo TDDFT – TDDFT module within Yambo, primarily for linear‑response optics. Community: CMP and optical properties. Ecosystem: post‑processing code reading DFT outputs (QE, etc.).
yambo-code
exciting TDDFT – linear‑response and real‑time TDDFT within exciting. Community: CMP, optical spectroscopies. Ecosystem: exciting FP‑LAPW code with GW and TDDFT modules.
iopscience.iop
+2
2.2 MBPT: GW / BSE / beyond
Yambo – MBPT code for GW and BSE; also TDDFT. Community: CMP and optical properties. Ecosystem: open source, interfaced to QE and other DFT codes.
yambo-code
BerkeleyGW – GW and GW‑BSE code for molecules and solids. Community: CMP, especially GW/BSE benchmarks. Ecosystem: open source, interfaced with many DFT codes.
berkeleygw
+1
ABINIT GW/BSE – ABINIT includes GW and Bethe–Salpeter implementations. Community: CMP and electronic excitations. Ecosystem: integrated in ABINIT workflow; used with Phonopy and post‑processors.
abinit
+1
exciting GW – GW implementation within exciting; includes G0W0 and beyond. Community: CMP and excited states. Ecosystem: part of exciting, which focuses on excitations.
exciting-code
+1
VASP GW – GW and (in newer versions) GW‑BSE features. Community: materials scientists with VASP licenses. Ecosystem: within VASP, integrated with PAW DFT and MBPT workflows.
vasp
+1
Spex – MBPT code (GW, GT, GWT, EELS, BSE) built on FLAPW basis (FLEUR family). Community: CMP, correlated materials. Ecosystem: part of MaX/FLAPW ecosystem; reads FLEUR electronic structure.
helmholtz
FHI-aims GW – all‑electron G0W0 implementation in FHI-aims. Community: CMP and quantum chemistry. Ecosystem: uses NAO basis and is integrated into FHI-aims workflows.
link.aps
Notes, cross‑check, uncertainties:

Pass 1: Wikipedia explicitly marks TDDFT and GW for some codes (e.g., Yambo, ABINIT, Octopus, exciting).
Pass 2: Additional TDDFT (SALMON, CP2K RT‑TDDFT) and more recent GW implementations (FHI‑aims GW, Spex, VASP GW) are added from official sites and literature.
salmon-tddft
+3
Pass 3: Many of these codes have explicit interfaces to ASE/AiiDA (QE, ABINIT, FHI-aims), or are part of the MaX ecosystem (FLEUR/Spex).
vasp
+2
Pass 4: EPW (electron–phonon) and Yambo/BerkeleyGW are key in excitonic/transport workflows; no conflicts.
Pass 5: Niche RT‑TDDFT research codes (e.g., ARTED, GCEED sub‑summed in SALMON) likely exist but are hard to enumerate exhaustively; completeness in RT‑TDDFT is uncertain.
3. Strongly correlated & many‑body methods
3.1 DMFT & beyond
Core DMFT toolkits & DFT+DMFT frameworks

TRIQS – Toolbox for Research on Interacting Quantum Systems (C++/Python library) with impurity solvers and DMFT components. Community: DMFT/strong correlations. Ecosystem: dft_tools for DFT+DMFT interfaces (Wien2k, VASP, etc.).
triqs.github
+1
TRIQS/cthyb – hybridization‑expansion CT‑QMC impurity solver (official TRIQS application). Community: DMFT impurity solvers. Ecosystem: used within TRIQS DMFT workflows and interfaced via dft_tools.
triqs.github
+1
DCore – integrated DMFT software for correlated electrons (DMFT workflow tools). Community: Japanese DMFT community and broader. Ecosystem: interfaces to DFT codes; uses TRIQS as backend.
TRIQS/dft_tools – DFT+DMFT interface between TRIQS and DFT codes (Wien2k, VASP, etc.). Community: DMFT practitioners. Ecosystem: AiiDA plugins; heavily used for Wien2k+DMFT.
github
Impurity solvers & DMFT codes

w2dynamics – hybridization‑expansion CT‑QMC impurity solver package (Wien/Würzburg). Community: DMFT impurity problems and DMFT+DFT. Ecosystem: works with DFT codes; standalone Fortran/C++ package.
github
+1
iQIST – quantum impurity solver toolkit (CT‑HYB and Hirsch‑Fye QMC). Community: DMFT impurity solvers. Ecosystem: open source; used in Asian and global DMFT community.
arxiv
+2
EDMFTF / eDMFT (Rutgers) – Haule’s DMFT implementation (eDMFT) interfaced to Wien2k and other DFT codes. Community: correlated materials. Ecosystem: embedded DMFT Functional; code and tutorials on Haule’s site.
hauleweb.rutgers
+2
ComDMFT – ab initio code combining LQSGW or DFT with DMFT, using Wannier90 and cRPA; part of Comscope. Community: correlated materials. Ecosystem: interfaces with Wannier90 and uses ComCTQMC solver.
github
+2
ComCTQMC – GPU‑accelerated CT‑QMC impurity solver (part of Comscope). Community: DMFT impurity solvers. Ecosystem: used as solver in ComDMFT workflows.
sciencedirect
+1
ComRISB – RISB/Gutzwiller component within ComDMFT. Community: strongly correlated systems, Gutzwiller approximations. Ecosystem: used together with ComDMFT for DFT+G/RISB.
github
DFT+DMFT implementations in DFT codes:
ABINIT – includes DMFT capabilities (charge self‑consistent DFT+DMFT).
abinit
Wien2k – DMFT via DMFTwWien2k and/or TRIQS/dft_tools; documented in Haule’s eDMFT tutorials.
hauleweb.rutgers
+1
Questaal – LMTO DFT code interfaced to Haule’s CTQMC for DFT+DMFT.
Other related many‑body toolkits

ALF (Algorithms for Lattice Fermions) – auxiliary‑field QMC code for lattice models (finite‑T and projective AFQMC). Community: lattice many‑body and model Hamiltonians. Ecosystem: open source; model‑focused but sometimes combined with DFT downfolded models.
gitpages.physik.uni-wuerzburg
+2
Notes, cross‑check, uncertainties:

Pass 1: DMFT is not covered systematically on Wikipedia; the baseline mentions only ABINIT and Yambo in DFT/MBPT context.
mattermodeling.stackexchange
Pass 2: Additional DMFT ecosystem (TRIQS, DCore, w2dynamics, iQIST, EDMFTF, ComDMFT, Questaal) is taken from official docs and Matter Modeling SE.
triqs.github
+3
Pass 3: TRIQS/dft_tools and AiiDA plugins explicitly provide DFT+DMFT workflows, particularly with Wien2k and VASP.
github
+1
Pass 4: DMFT workflows heavily rely on Wannier90 and DFT codes; this is consistent with tools in Sections 1 and 5.
Pass 5: There are many smaller DMFT and impurity solver codes, especially model‑only codes and experimental ones; completeness for DMFT solvers cannot be guaranteed.
3.2 Quantum Monte Carlo (QMC)
QMCPACK – high‑performance QMC code for VMC, DMC, AFQMC, reptation QMC. Community: QMC in CMP and quantum chemistry. Ecosystem: open source, HPC‑oriented, interfaces to Quantum ESPRESSO, NWChem, PySCF, etc.
qmcpack
+2
CASINO – QMC code for VMC and DMC. Community: QMC community (originally Cambridge). Ecosystem: academic license; widely used for solids and molecules.
en.wikipedia
+1
TurboRVB – VMC and lattice‑regularized DMC (LRDMC) with RVB wavefunctions. Community: strongly correlated systems, superconductors. Ecosystem: open source, supports geometry optimization and forces for MD.
github
+2
ALF – AFQMC for lattice fermions (finite‑T and projective). Community: lattice models and strongly correlated electrons. Ecosystem: open source, model‑Hamiltonian oriented.
gitpages.physik.uni-wuerzburg
+1
(UNCERTAIN) There are other AFQMC and phaseless AFQMC research codes (e.g., group‑specific codes) not cataloged centrally; completeness in QMC cannot be guaranteed.
Notes, cross‑check:

Pass 1: CASINO and QMCPACK appear on Wikipedia.
Pass 2: TurboRVB and ALF added from official sources.
github
+1
Pass 3: QMCPACK interfaces with QE/PySCF and is used in QMCPACK summer schools; CASINO and TurboRVB have their own ecosystems.
qmcpack
+1
Pass 4: QMC codes provide input data for transport and topology studies (e.g., excited‑state QMC); no conflicts.
Pass 5: Specialized QMC codes (e.g., for specific lattice models) are numerous; not all are listed.
4. Wavefunction‑based quantum chemistry (many‑body)
4.1 Coupled cluster methods
ORCA – general‑purpose quantum chemistry with emphasis on spectroscopy; includes CC methods. Community: quantum chemistry, spectroscopists. Ecosystem: academic license; widely used in open‑shell chemistry.
faccts
+1
PSI4 – open‑source CC/DFT code, highly modular. Community: quantum chemistry methods development. Ecosystem: open source (LGPL), Python interface, integrated with various workflows.
psicode
+1
PySCF – Python quantum chemistry library with CC, MP2, CASSCF, etc. Community: quantum chemistry and methods. Ecosystem: open source (Apache), used in high‑throughput and ML workflows.
pyscf
+1
CFOUR – Coupled‑Cluster techniques for Computational Chemistry. Community: high‑accuracy CC benchmarks. Ecosystem: academic license, focused on post‑HF and CC.
cfour.uni-mainz
MRCC – suite of CC/CI programs, including multi‑reference CC. Community: high‑accuracy quantum chemistry. Ecosystem: academic/commercial; interfaces to Cfour, Columbus, Dirac, Molpro, ORCA, PSI.
mrcc
+2
Molpro – ab initio suite emphasizing CC, MRCI, and explicitly correlated methods. Community: high‑accuracy quantum chemistry. Ecosystem: commercial license; extensive CC and local‑CC methods.
molpro
+1
Dalton – suite for molecular properties with HF, DFT, MCSCF, CC. Community: molecular properties and response. Ecosystem: open source (LGPL); part of Dalton/LSDalton ecosystem.
daltonprogram
+2
DIRAC – relativistic ab initio code with HF, DFT, CC. Community: relativistic quantum chemistry. Ecosystem: open source, specialized in 4‑component and 2‑component relativistic methods.
diracprogram
+1
4.2 Configuration interaction & multireference
OpenMolcas – multiconfigurational quantum chemistry package (CASSCF, CASPT2, etc.). Community: multireference quantum chemistry. Ecosystem: open source; supports pair‑density functional theory and dynamics.
molcas.gitlab
+1
BAGEL – parallel quantum chemistry package with multireference methods (CASPT2 gradients, relativistic MR, etc.). Community: multireference and excited‑state dynamics. Ecosystem: open source; uses SMITH3 code generator.
arxiv
+1
COLUMBUS – multireference CI program (MR‑CISD, MR‑AQCC) with analytic gradients. Community: multireference dynamics and nonadiabatic processes. Ecosystem: open source, used in photochemistry.
columbus-program-system.gitlab
+2
PySCF MR – multireference capabilities (CASSCF, NEVPT2, etc.). Community: Python‑based MR community. Ecosystem: integrated with PySCF’s Python ecosystem.
en.wikipedia
Notes, cross‑check:

Pass 1: Wikipedia includes Dalton, DIRAC, ORCA, PSI4, PySCF, Molpro, OpenMolcas, MRCC, COLUMBUS, etc., in its tables.
Pass 2: BAGEL and newer PySCF developments are from their official pages/papers; DIRAC and Dalton are from their official documentation.
arxiv
+1
Pass 3: PSI4, PySCF, pymatgen, ASE, and workflow engines (AiiDA, atomate) are well integrated; Dalton/Dirac have their own tooling.
Pass 4/5: Wavefunction codes mainly serve molecular clusters and molecules rather than periodic solids, though PySCF has periodic capabilities. Completeness for specialized MR codes is uncertain (many small, group‑specific codes exist).
5. Tight‑binding, model Hamiltonians & downfolding
Wannier90 – code for maximally localized Wannier functions and Wannier interpolations. Community: CMP/MS. Ecosystem: heavily used across codes (VASP, QE, FLEUR, etc.) for transport, topology, e‑ph.
wannier
+1
WannierTools – topological materials package operating on TB models from Wannier90. Community: topological materials. Ecosystem: reads Wannier90 HR files; computes surface states, Wilson loops, Weyl/Dirac points.
wanniertools
+1
PythTB – Python tight‑binding package for constructing and analyzing TB models (band structures, topology). Community: topology and model‑Hamiltonian studies. Ecosystem: open source, educational and research use.
pythtb.readthedocs
+2
TBmodels – Python package for evaluating TB models and symmetrizing Wannier‑like models. Community: transport and topology. Ecosystem: often used with Wannier90; integrated with AiiDA tbextraction workflows.
tbmodels.greschd
TRIQS/dft_tools – DFT+DMFT interface that generates effective lattice models (impurity and lattice) from DFT via Wannier downfolding. Community: DMFT model building. Ecosystem: TRIQS/AiiDA ecosystem.
github
BoltzWann – module/part of Wannier90 ecosystem for Boltzmann transport using Wannier functions. Community: thermoelectrics and transport. Ecosystem: tightly integrated with Wannier90.
wannier90.readthedocs
+2
Notes, cross‑check:

Pass 1: Some TB tools appear on Wikipedia (Wannier90, WannierTools; PythTB not as prominently).
wannier
Pass 2: PythTB, TBmodels, BoltzWann, and TRIQS/dft_tools added from official docs and community pages.
pythtb.readthedocs
+2
Pass 3: Wannier90 and TBmodels have AiiDA plugins (aiida‑tbextraction); PythTB often used with ASE.
Pass 4: Phonon/transport/topology tools heavily rely on Wannier90 and TBmodels; EPW, Z2Pack, WannierTools, etc., are consistent.
Pass 5: Numerous small TB codes for specific systems (2D materials, specific materials classes) exist; completeness for TB tools is not guaranteed.
6. Phonons, lattice dynamics & electron–phonon
Harmonic / quasi‑harmonic phonons

Phonopy – harmonic and quasi‑harmonic phonon calculations (supcell finite‑displacements). Community: CMP/MS. Ecosystem: interfaces to VASP, QE, CRYSTAL, ABINIT, DFTB+, etc.
phonopy.github
+2
ABINIT phonons (DFPT) – ABINIT includes DFPT for phonons, IR/Raman, Born effective charges. Community: CMP. Ecosystem: part of ABINIT; workflows can be automated via ASE/AiiDA.
docs.abinit
+1
Quantum ESPRESSO PHonon – PHonon module in QE for DFPT phonons. Community: QE users. Ecosystem: integrated in QE; often combined with Phonopy for supercell calculations.
phonopy.github
Anharmonic phonons & lattice thermal conductivity

Phono3py – phonon–phonon interaction and lattice thermal conductivity (RTA and direct BTE). Community: CMP thermal transport. Ecosystem: interfaces to VASP, QE, CRYSTAL, TURBOMOLE, ABINIT.
phonopy.github
+1
ALAMODE – anharmonic lattice dynamics and thermal conductivity code. Community: lattice dynamics and thermal transport. Ecosystem: uses DFT forces from VASP, QE, OpenMX, LAMMPS, etc.
alamode.readthedocs
+1
TDEP – Temperature Dependent Effective Potential code for finite‑T lattice dynamics. Community: anharmonic phonons at finite temperature. Ecosystem: Fortran/C++ package, integrates with DFT force calculations.
joss.theoj
+3
ShengBTE – solver of phonon Boltzmann transport equation for lattice thermal conductivity. Community: thermal transport. Ecosystem: open source, uses second/third‑order force constants.
shengbte
+2
almaBTE – alternative BTE solver for lattice thermal conductivity. Community: thermal transport. Ecosystem: complementary to ShengBTE, used in similar workflows.
shengbte
+1
THERMACOND / ALATDYN – codes for lattice thermal conductivity and thermodynamic properties based on anharmonic force constants. Community: thermal transport and anharmonicity. Ecosystem: newer research codes; THERMACOND for phonon transport.
sciencedirect
+1
SSCHA – Stochastic Self‑Consistent Harmonic Approximation for strongly anharmonic phonons. Community: strongly anharmonic systems and phase diagrams. Ecosystem: Python package, interfaces to DFT codes.
sscha
+2
Electron–phonon coupling & superconductivity

EPW – Electron–Phonon Wannier code for e‑ph coupling, transport, superconductivity (Migdal–Eliashberg). Community: e‑ph and superconductors. Ecosystem: uses DFPT + Wannier90; interfaced with QE and ABINIT.
docs.epw-code
+2
ABINIT e-ph – ABINIT supports DFPT‑based electron–phonon coupling calculations. Community: CMP and superconductivity. Ecosystem: integrated in ABINIT workflows.
abinit
+1
Notes, cross‑check:

Pass 1: Wikipedia mentions some phonon tools (e.g., Phonopy), but not most; ShengBTE, Phono3py appear only in specialist contexts.
Pass 2: Added Phono3py, ALAMODE, TDEP, ShengBTE, almaBTE, SSCHA, EPW, ABINIT e-ph, and THERMACOND/ALATDYN from their official sites and literature.
phonopy.github
+3
Pass 3: Phonopy, phono3py explicitly support VASP, QE, CRYSTAL, ABINIT, DFTB+; TDEP and SSCHA interfaced to DFT codes; EPW uses DFPT + Wannier90.
phonopy.github
+2
Pass 4: TDEP, ALAMODE, ShengBTE, almaBTE, EPW, SSCHA form coherent thermal‑transport/e‑ph ecosystem.
Pass 5: Smaller research codes for anharmonic phonons or specific materials likely exist; completeness here is uncertain.
7. Molecular & ab initio dynamics
CPMD – Car–Parrinello MD code based on plane‑wave DFT. Community: liquids, interfaces, CMP. Ecosystem: open source; interfaces to i‑PI for advanced sampling.
cpmd
+1
CP2K MD – BO/MD, CPMD, metadynamics, path‑integral MD (PINT). Community: condensed matter and chemistry. Ecosystem: PINT module and i‑PI interface support advanced PIMD.
cp2k
+2
VASP MD – ab initio MD using VASP’s DFT engine. Community: materials science. Ecosystem: within VASP; used in high‑throughput workflows (Materials Project).
vasp
+1
Quantum ESPRESSO MD – CP (Car–Parrinello) and Verlet MD modules. Community: CMP. Ecosystem: within QE; used in Phonopy and EPW workflows.
phonopy.github
i-PI – universal Python interface for path‑integral MD and advanced sampling, acting as server to DFT codes. Community: nuclear quantum effects in solids and molecules. Ecosystem: interfaces to CP2K, FHI-aims, DFTB+, etc.
gle4md
+3
ABINIT PIMD – ABINIT includes path‑integral MD capabilities. Community: quantum nuclei effects. Ecosystem: integrated in ABINIT.
abinit
Notes, cross‑check:

Pass 1: CPMD and CP2K MD are on Wikipedia; i-PI is not.
Pass 2: Added i‑PI as key MD wrapper; ABINIT PIMD from ABINIT docs.
gle4md
+1
Pass 3: i-PI is widely interfaced to CP2K, FHI-aims, DFTB+, etc., forming a cross‑code MD ecosystem.
dftbplus-recipes.readthedocs
Pass 4: Phonon codes often rely on MD or force‑constant data from these MD engines; EPW and TDEP explicitly use DFT outputs.
Pass 5: Many advanced MD methods (ring‑polymer MD, advanced PIMD variants) have smaller group codes; completeness for MD cannot be guaranteed.
8. Structure prediction & global optimization
USPEX – evolutionary algorithm for crystal structure prediction (Universal Structure Predictor: Evolutionary Xtallography). Community: crystal structure prediction. Ecosystem: standalone, used for 0–3D systems; interfaced to DFT codes.
uspex-team
+1
XtalOpt – open‑source multi‑objective evolutionary crystal structure prediction code. Community: materials discovery. Ecosystem: interfaces to VASP, QE, ABINIT, etc.; BSD license.
xtalopt.github
+1
CALYPSO – particle‑swarm optimization method (Crystal Structure AnaLYsis by Particle Swarm Optimization). Community: structure prediction. Ecosystem: standalone, interfaced to external DFT codes.
calypso
+1
AIRSS – ab initio random structure searching. Community: structure prediction. Ecosystem: simple scripts, highly parallel, often used with ABINIT and other DFT codes.
airss-docs.github
+1
GASP – genetic algorithm for structure prediction (part of the XtalOpt family in earlier works). Community: structure prediction and evolutionary searches. Ecosystem: integrated into XtalOpt ecosystem and used with external DFT codes.
xtalopt.github
Notes, cross‑check:

Pass 1: USPEX, XtalOpt, CALYPSO, AIRSS are not systematically covered by Wikipedia’s baseline; they are specialist codes known within the crystal‑prediction community.
Pass 2: Added descriptions from their official pages and key papers.
uspex-team
+3
Pass 3: AIRSS and XtalOpt are often combined with VASP/QE/ABINIT; XtalOpt documented as interfaced to multiple codes; AIRSS documented as “ab initio random structure searching”.
airss-docs.github
+1
Pass 4: These structure prediction tools heavily depend on DFT and phonon codes; no conflicts.
Pass 5: Additional local structure prediction or surface search tools (e.g., for surfaces/clusters) are numerous; completeness for structure prediction cannot be guaranteed.
9. Post‑processing, analysis & visualization
Electronic structure, DOS, bands, topology

PyProcar – Python library to plot band structures and Fermi surfaces from PROCAR files; supports VASP, Elk, QE, ABINIT, LOBSTER. Community: CMP/MS. Ecosystem: open source; used for band unfolding and projections.
pyprocar.readthedocs
+1
FermiSurfer – Fermi‑surface visualization tool. Community: CMP/metals. Ecosystem: reads outputs from multiple DFT codes.
mitsuaki1987.github
+2
WannierTools – topological analysis (Wilson loops, surface states) from TB models. Community: topological materials. Ecosystem: uses Wannier90 outputs; widely used for topology.
wanniertools
+1
Z2Pack – code for computing topological invariants (Z2, Chern numbers, etc.). Community: topological insulators/semimetals. Community: topology; ecosystem: works with TB and DFT Hamiltonians via Wannier90.
docs.abinit
+2
Sumo – Python toolkit for plotting and analysis of ab initio solid‑state data (band structures, DOS, optical absorption). Community: CMP visualization. Ecosystem: open source, builds on spglib and pymatgen/ASE.
github
Transport and thermoelectric properties

BoltzTraP – semiclassical Boltzmann transport properties from band structures. Community: thermoelectrics and transport. Ecosystem: interfaced to WIEN2k, ABINIT, SIESTA, VASP, QE.
wiki.rc.usf
+1
BoltzWann – Boltzmann transport using Wannier functions. Community: thermoelectrics. Ecosystem: part of Wannier90 ecosystem.
wannier90.readthedocs
+1
Bonding and chemical analysis

LOBSTER – Local Orbital Projections, Atomic Charges, and Chemical Bonding Analysis from PAW‑based DFT. Community: CMP and quantum chemistry. Ecosystem: supports VASP and other PAW codes; LobsterPy provides automated post‑processing.
psi-k
+3
LobsterPy – Python package to automate LOBSTER analysis and generate ML features. Community: CMP/ML. Ecosystem: depends on LOBSTER outputs.
github
DFT‑specific post‑processing kits

Vaspkit – pre‑ and post‑processing program for VASP (k‑paths, symmetry, elastic, band structures, etc.). Community: VASP users. Ecosystem: standalone, focused on VASP but extensible.
vaspkit
+3
Notes, cross‑check:

Pass 1: Many of these tools (PyProcar, FermiSurfer, BoltzTraP, Vaspkit, Sumo) are not on the Wikipedia baseline; they are community codes often cited in papers.
pyprocar.readthedocs
+2
Pass 2: Added PyProcar, FermiSurfer, Sumo, Lobster/LobsterPy, Vaspkit, Z2Pack, BoltzWann from their official docs and papers.
pyprocar.readthedocs
+3
Pass 3: PyProcar, Sumo, BoltzTraP, and Vaspkit build on ASE/pymatgen or support their I/O; Z2Pack and WannierTools interface Wannier90; AiiDA supports some analysis workflows.
Pass 4: These tools are designed specifically for phonon/transport/topology workflows and are consistent with DFT/MBPT outputs.
Pass 5: Many smaller analysis utilities exist (e.g., specialized visualization tools for certain codes). Completeness in analysis/visualization cannot be guaranteed.
10. Frameworks, workflow engines & databases
Core Python libraries and workflow engines

ASE (Atomic Simulation Environment) – Python library for setting up, running, and analyzing atomistic simulations; includes calculator interfaces to many DFT/MD/TB codes. Community: general CMP/MS. Ecosystem: foundational library used across the ecosystem.
ase-lib
+1
pymatgen – Python Materials Genomics library (structures, analysis, I/O, defect/phase diagram workflows). Community: materials informatics. Ecosystem: integrates with ASE, Materials Project, atomate, AiiDA.
pymatgen
+1
FireWorks – workflow system for defining and executing scientific workflows (Python/JSON/YAML, MongoDB backend). Community: high‑throughput computational science. Ecosystem: basis for atomate; used in Materials Project.
legacy.materialsproject
+2
custodian – job‑management and error‑correction library (especially for high‑throughput calculations). Community: high‑throughput materials. Ecosystem: used within atomate and Materials Project workflows.
materialsproject.github
High‑level workflow layers

atomate – high‑level materials workflows built on pymatgen, FireWorks, and custodian (VASP‑centric). Community: high‑throughput DFT workflows. Ecosystem: open source; predecessor to atomate2.
sciencedirect
+1
atomate2 – modular workflow framework for materials science, multi‑code support (VASP, QE, etc.). Community: high‑throughput CMP. Ecosystem: extends atomate design; integrated with custodian and pymatgen.
pmc.ncbi.nlm.nih
Provenance and data infrastructure

AiiDA – Automated Interactive Infrastructure and Database for Computational Science. Community: CMP/MS and computational science. Ecosystem: plugin registry with >100 plugins for many codes (QE, VASP, CP2K, Wannier90, Yambo, etc.).
aiida
+3
AiiDA plugin registry – public list of AiiDA plugins (aiida‑quantumespresso, aiida‑vasp, aiida‑cp2k, aiida‑wannier90, etc.). Community: AiiDA users. Ecosystem: official registry at aiida.net/plugins.
aiida
+1
Materials Project – high‑throughput DFT database and associated tools (pymatgen, atomate). Community: materials informatics. Ecosystem: depends on pymatgen, atomate, custodian, FireWorks; integrates with AiiDA workflows.
pymatgen
+1
Materials Cloud – repository and platform powered by AiiDA, storing data and workflows. Community: European and global CMP. Ecosystem: built on AiiDA; hosts AiiDA workflows and AiiDAlab apps.
materialscloud
+2
Additional infrastructure

libAtoms/Quippy – extension of ASE for more advanced potentials and interfaces (used in some chemistry and HPC contexts). Community: ASE users. Ecosystem: part of the ASE ecosystem (not detailed above due to size).
ase-workshop-2023.github
MDI / MolSSI drivers – cross‑code driver standard (e.g., MDI interface between LAMMPS and QM codes). Community: computational chemistry and materials. Ecosystem: used to connect LAMMPS with VASP, CP2K, etc.
docs.lammps
+1
Notes, cross‑check:

Pass 1: Wikipedia lists neither ASE nor pymatgen/FireWorks/AiiDA explicitly; it is limited to “simulation codes”, not infrastructures.
Pass 2: ASE, pymatgen, FireWorks, custodian, atomate, atomate2, AiiDA, Materials Cloud added from official docs and key papers.
ase-lib
+4
Pass 3: ASE and pymatgen are foundational; atomate/atomate2 build on FireWorks/custodian/pymatgen; AiiDA provides its own plugin ecosystem; Materials Cloud and Materials Project are major databases.
Pass 4: All these frameworks are designed to orchestrate DFT, phonon, MBPT, DMFT, QMC, etc., and are consistent with the tools listed in Sections 1–9.
Pass 5: Many smaller workflow managers and lab‑specific pipelines exist; completeness for workflow engines is difficult to guarantee.
11. Small, niche, community & research‑grade tools (selected examples)
Per your instruction, I cannot guarantee that every niche tool is included, but here are selected tools that are “small”, research‑grade, or frequently omitted from generic lists, with verified sources:

DFT‑F – adaptive finite‑element DFT for large‑scale calculations (HPC).
sciencedirect
NESSIE – FEM‑based all‑electron DFT + RT‑TDDFT.
nessie-code
DMFTwDFT – DFT+DMFT project (not covered in detail above; appears in the DMFT SE discussion). Community: DMFT. Ecosystem: likely integrates DFT codes with DMFT solvers (UNCERTAIN beyond this).
ALF – AFQMC for lattice fermions (already in QMC).
scipost
w2dynamics – CT‑HYB impurity solver (already in DMFT).
github
ComCTQMC – GPU‑accelerated CT‑QMC impurity solver (part of Comscope).
sciencedirect
+1
SALMON – RT‑TDDFT for molecules and solids (already in TDDFT).
salmon-tddft
BoltzWann – Wannier‑based Boltzmann transport (already in TB).
sciencedirect
Simphony – Wannier‑based package for topological phonons (phonon topology). Community: phonon topology. Ecosystem: research‑grade but publicly documented; extends Wannier ecosystem to phonons.
arxiv
Z2Pack – topological invariants from DFT/TB models.
docs.abinit
+1
WannierTools – topological materials package (also in TB).
wanniertools
TBmodels – Wannier‑like TB model evaluation/symmetrization.
tbmodels.greschd
Sumo – solid‑state visualization toolkit.
github
LobsterPy – automation of LOBSTER bonding analysis.
github
Uncertainties and gaps:

Many research codes are single‑group or single‑project tools and lack stable websites or DOIs; I cannot reasonably enumerate them all.
Some tools listed on the PSI‑k and NOMAD code lists are not captured above in detail (e.g., some more specialized codes on psi‑k “software” page).
psi-k


Computational Tools in Condensed Matter Physics, Materials Science, Quantum Chemistry, and Solid‑State Physics
Purpose: A maximally complete, structured enumeration of computational tools used in:

Condensed matter physics
Materials science
Quantum chemistry
Solid‑state physics
Note: Completeness cannot be mathematically guaranteed. This list draws from multipleauthoritative sources (Wikipedia tables, official code documentation, community lists, peer‑reviewedpapers). I have not guessed tool names and have marked uncertain assignments where appropriate.

Method classes
Ground‑State Electronic Structure (DFT & Variants)
Time‑Dependent & Excited‑State Methods
Strongly Correlated & Many‑Body Methods
Wavefunction‑Based Quantum Chemistry (Many‑Body)
Tight‑Binding, Model Hamiltonians & Downfolding
Phonons, Lattice Dynamics & Electron–Phonon
Molecular & Ab Initio Dynamics
Structure Prediction & Global Optimization
Post‑Processing, Analysis & Visualization
Frameworks, Workflow Engines & Databases
Small, Niche, Community & Research‑Grade Tools
1. Ground‑State Electronic Structure (DFT & Variants)
1.1 Plane‑wave / real‑space / wavelet / mixed (mainly solids)
Quantum ESPRESSO – Plane‑wave pseudopotential/PAW DFT, hybrids, DFPT, MD.
VASP – Plane‑wave PAW DFT, hybrids, meta‑GGA, GW/BSE, RPA/ACFDT, MD.
ABINIT – Plane‑wave pseudopotential DFT, TDDFT, GW/BSE, phonons/e‑ph (DFPT, EPH).
CASTEP – Plane‑wave pseudopotential DFT for solids (academic license).
CP2K – Mixed Gaussian + plane‑waves (GPW/GAPW) DFT, MP2, RPA, AIMD, metadynamics.
GPAW – Real‑space multigrid PAW DFT; supports TDDFT and GW.
BigDFT – Wavelet‑basis DFT.
Qbox – Plane‑wave pseudopotential DFT.
PEtot – Plane‑wave pseudopotential DFT for large systems.
DACAPO – Plane‑wave pseudopotential DFT.
Socorro – Plane‑wave pseudopotential DFT.
PARSEC – Real‑space pseudopotential DFT.
SPHINX – Plane‑wave DFT code (Jülich family).
ONETEP – Linear‑scaling DFT with NGWFs (periodic solids).
OpenAtom – Plane‑wave DFT with Charm++ parallelization.
CONQUEST – Linear‑scaling DFT with NAO/spline bases, periodic solids.
RESCU – Grid/NAO/PW electronic structure for solids.
1.2 Localized / Gaussian / numeric atom‑centered basis (solids + molecules)
SIESTA – Numerical atom‑centered orbital (NAO) basis pseudopotential DFT for solids.
OpenMX – NAO basis pseudopotential DFT (open source).
FHI-aims – All‑electron numeric atom‑centered orbital basis, solids and molecules; hybrids, GW/BSE.
CRYSTAL – Gaussian‑type orbital all‑electron code for 0D–3D periodic systems; HF, DFT, hybrids.
BAND (part of AMS/ADF) – LCAO periodic DFT with Slater/numerical orbitals; supports TDDFT for solids.
ADF (Amsterdam Density Functional) – Molecular DFT with Slater‑type orbitals; part of the Amsterdam Modeling Suite.
DMol3 – Numerical atomic‑orbital DFT (Materials Studio).
ORCA – GTO basis molecular electronic structure; strong in correlated wavefunction methods; also supports DFT and hybrids.
Gaussian – GTO molecular QC; HF, DFT (including hybrids), post‑HF methods.
Q‑Chem – Molecular DFT and wavefunction methods (HF, CC, TDDFT).
NWChem – Molecular and periodic (PW) electronic structure; HF, DFT, CC, TDDFT (real‑time in some modules).
PSI (PSI4) – Molecular QC; HF, DFT, CC, TDDFT.
PySCF – Python‑based molecular QC; HF, DFT, CC, multireference, GW/BSE modules.
Dalton – HF and response properties; MP2, CC, etc.
DIRAC – Relativistic molecular QC (HF, DFT, CC, etc.).
Molpro – Molecular HF, DFT, CC, multireference methods.
OpenMolcas / MOLCAS – Multireference CASSCF/CASPT2 with DFT and post‑HF capabilities.
MRCC – Molecular CC and multireference methods (CCSD, CCSD(T), etc.).
MPQC – Massively Parallel Quantum Chemistry (HF, DFT, CC, etc.).
ACES – Coupled‑cluster code; primarily molecular.
COLUMBUS – Molecular multireference CI and dynamics.
eT – Molecular electronic structure (HF, CC).
NESSIE – Finite‑element DFT code.
MADNESS – Multiresolution adaptive numerical environment for electronic structure; real‑space wavelets.
Octopus – Real‑space grid code; focused on TDDFT and electron dynamics, but also supports ground‑state DFT.
PARSEC – Real‑space pseudopotential DFT.
1.3 Augmented methods / all‑electron periodic
WIEN2k – Full‑potential LAPW code for solids.
FLEUR – FLAPW code (Jülich family; open source).
ELK – FP‑LAPW code (GPL).
exciting – FP‑LAPW with TDDFT, MBPT, GW/BSE and phonon/e‑ph modules.
RSPt – FP‑LMTO code.
FPLO – Full‑potential local‑orbital code.
AIMPRO – Local‑orbital DFT code (for defects).
1.4 Semi‑empirical / tight‑binding DFT (ground‑state approximate)
DFTB+ – Density Functional Tight Binding implementation for molecules and solids; fast approximate DFT.
DFTB module in AMS (ADF/BAND) – Integration of DFTB into Amsterdam Modeling Suite workflow.
xtb – Extended tight‑binding (GFN‑xTB etc.) package for fast approximate electronic structure.
1.5 Hybrids / exact exchange in ground‑state DFT
VASP – Implements many hybrid functionals (HSE03, HSE06, PBE0, etc.) and ACFDT/RPA.
Quantum ESPRESSO – Supports local, semi‑local, and hybrid DFT functionals; improved parallel treatment of exact exchange.
CRYSTAL – All‑electron periodic code with hybrid functionals (global and range‑separated).
FHI-aims – All‑electron NAO code with highly optimized hybrid DFT capabilities.
ORCA, Q‑Chem, Gaussian, Molpro, Dalton, DIRAC, PSI4, PySCF – Provide hybrid functionals for molecular systems (well‑known in quantum chemistry).
2. Time‑Dependent & Excited‑State Methods
2.1 TDDFT & real‑time propagation
Octopus – Real‑space real‑time TDDFT for molecules and solids; also supports Ehrenfest dynamics.
GPAW – Real‑time TDDFT module (PAW).
NWChem – Has real‑time TDDFT module (RTTDDFT).
ABINIT – Supports TDDFT (Casida formalism) for molecules.
Quantum ESPRESSO – TDDFT module (TD‑DeltaSCF / linear‑response) via turboTDDFT/TDDF modules.
FHI-aims – TDDFT capabilities; see tutorials on GW/BSE and TDDFT.
exciting – Module for time‑dependent DFT as part of its excited‑state focus.
eT – Molecular CC with time‑dependent capabilities (TDHF).
MOLCAS / OpenMolcas, Dalton, DIRAC – Provide linear‑response TDHF/TDDFT in molecular codes.
2.2 MBPT (GW / BSE / beyond)
Yambo – GW and BSE code interfaced mainly with Quantum ESPRESSO.
BerkeleyGW – GW and GW+BSE code for molecules and solids.
Spex – MBPT code based on FLAPW; GW, COHSEX, BSE, EELS, RPA total energy; part of the Jülich FLEUR family.
West – GPU‑accelerated code for full‑frequency GW and related properties.
VASP – GW and BSE modules, plus RPA/ACFDT.
ABINIT – GW implementation; G₀W₀, partially self‑consistent variants, and BSE support (via Bethe–Salpeter notes).
exciting – Focus on excited states; supports GW and BSE.
CP2K – GW implementation; CP2K tutorial describes GW for molecules and solids.
FHI-aims – GW implementation and GW+BSE for molecules/solids.
Fiesta – Gaussian‑basis GW and BSE code; reads states from ORCA/NWChem.
Elk – FP‑LAPW code with advanced features; has been used in GW contexts (UNCERTAIN if GW is a core feature).
2.3 RPA / beyond‑DFT correlation (ground‑state but beyond standard DFT)
VASP – Implements RPA and ACFDT correlation energy (ALGO=RPA).
Spex – Supports RPA total‑energy calculations.
CP2K – RPA via RI‑RPA/ACFDT modules.
ABINIT – RPA‑related RPACorrEn topic; supports correlation‑energy RPA‑style calculations.
Quantum ESPRESSO – Provides ACFDT/RPA via the qe‑acfdt module.
3. Strongly Correlated & Many‑Body Methods
3.1 DMFT & beyond
TRIQS – Toolbox for DMFT and many‑body solvers; core of TRIQS/DFTTools.
TRIQS/DFTTools – DFT+DMFT framework interfaced to Wien2k, Elk, VASP, and others for realistic materials.
DCore – Integrated DMFT software; built on TRIQS, ALPS; interfaces to Quantum ESPRESSO and OpenMX.
EDMFTF – DFT+DMFT code integrated with WIEN2k; used with impurity solvers like CT‑HYB.
eDMFT (DFT+Embedded DMFT Functional, Haule) – DFT+DMFT implementation using WIEN2k with CTQMC/OCA/NCA solvers.
ComDMFT – Ab initio GW+EDMFT code (v2.0 full GW+EDMFT). Built on Wannier90, ComCTQMC and LqsgwFlapw.
DMFTwDFT – DFT+DMFT code using the EDMFTF impurity solver.
w2dynamics – Continuous‑time QMC impurity solver for DMFT (CT‑HYB). Used in Wien2k+DMFT workflows.
iQIST – Open‑source CT‑QMC impurity solver toolkit for DMFT.
ALPS DMFT application – Provides CT‑QMC impurity solvers and DMFT self‑consistency loops.
ALPS core libraries – Generic libraries for strongly correlated lattice models and DMFT solvers.
AMULET – DFT+DMFT code interfaced to Quantum ESPRESSO and ELK.
Portobello – Framework comparing DFT+DMFT and LDA+G (Gutzwiller), used in comparative studies.
CyGutz, ComGutz – Open‑source LDA+RISB/Gutzwiller implementations.
3.2 Quantum Monte Carlo
3.2.1 Continuum QMC (VMC, DMC, AFQMC, etc.)
QMCPACK – Ab initio QMC for atoms, molecules, and solids; VMC, DMC, reptation Monte Carlo, and orbital‑space AFQMC.
CASINO – Continuum QMC code for molecules and periodic solids; VMC and DMC.
TurboRVB – QMC package; VMC and lattice‑regularized DMC (LRDMC) for molecules and solids.
CHAMP – Cornell–Leiden ab initio materials package; VMC, DMC, wave‑function optimization; supports molecules and solids.
QWalk – Quantum Monte Carlo program for molecules and solids (VMC, DMC, reptation).
3.2.2 Model Hamiltonian QMC / lattice QMC
ALPS / ALPS DMFT – CT‑QMC impurity solvers and lattice model QMC (Hubbard models, etc.).
CT‑INT / CT‑HYB (ALPS) – CT‑QMC impurity solvers for interactions and hybridization expansions.
ComCTQMC – CT‑QMC impurity solver used in ComDMFT.
4. Wavefunction‑Based Quantum Chemistry (Many‑Body)
4.1 Coupled cluster methods
ORCA – CCSD, CCSD(T), local CC, etc.
MRCC – Arbitrary‑order CC methods and multireference CC; also coupled cluster with singles, doubles (CCSD, CCSD(T)).
PSI / PSI4 – CCSD, CCSD(T), higher excitations in some cases, EOM‑CC.
PySCF – Python‑based CCSD, CCSD(T), and EOM‑CC; also FCI and multireference.
Molpro – CCSD, CCSD(T), and higher CC methods; MR‑CC for multireference.
Dalton – CC methods including CCSD and CCSD(T).
DIRAC – Relativistic CC and response methods.
ACES – High‑accuracy CC methods.
COLUMBUS – Focus on multireference CI; some CC capabilities through CI expansions.
Q‑Chem, Gaussian, CFOUR – Provide standard CC methods; CFOUR is a specialized high‑accuracy CC code.
4.2 Configuration interaction & multireference
OpenMolcas / MOLCAS – CASSCF, CASPT2, NEVPT2 and related multireference methods.
BAGEL – Multireference quantum chemistry (CASSCF, RAS, etc.).
PySCF – Supports CASSCF, DMRG, and related multireference methods.
Molpro – MRCI and related methods; multireference CI approaches.
COLUMBUS – Strong focus on CI, including MRCI; supports excited states.
MRCC – Multireference CI and coupled cluster (MR‑CC).
Dalton – Offers some CI and multiconfigurational capabilities.
5. Tight‑Binding, Model Hamiltonians & Downfolding
Wannier90 – Maximally localized Wannier functions (MLWFs); interfaced to many DFT codes.
WannierTools – TB and topological analysis code using Wannier90 TB models; Wilson loops, Weyl/Dirac points, surface states.
PythTB – Python tight‑binding package; models and topology (Berry phases, Wilson loops, etc.).
TBmodels – Python package for tight‑binding models; symmetrization and format conversion.
TRIQS/DFTTools – Constructs TB Hamiltonians via Wannier90 and implements DMFT self‑consistency.
DFTB+ – Tight‑binding DFT; constructs TB models from parameters.
xtb – Extended TB (GFN‑xTB) model Hamiltonian package.
Tight‑Binding Studio – Tool to construct TB models via LCAO and two‑center approximation from DFT.
6. Phonons, Lattice Dynamics & Electron–Phonon
Phonopy – Harmonic and quasi‑harmonic phonons; phonon band structures, DOS, thermal properties.
phono3py – Phonon–phonon interaction and lattice thermal conductivity via BTE; interfaces to VASP, QE, etc.
ShengBTE – Solver of the BTE for phonons; computes lattice thermal conductivity.
almaBTE – BTE solver for phonons; related to ShengBTE; thermal conductivity vs temperature.
THERMACOND – Code computing lattice thermal conductivity using phonon BTE; benchmarked vs ShengBTE.
TDEP – Temperature dependent effective potential method; builds force constants and phonons from AIMD trajectories.
ALAMODE – Anharmonic lattice dynamics; harmonic/anharmonic force constants and lattice thermal conductivity.
EPW – Electron–phonon coupling using DFPT and Wannier interpolation; e–ph linewidths, transport, superconductivity.
Quantum ESPRESSO phonons (PHonon) – DFPT lattice dynamics and electron–phonon coupling.
ABINIT DFPT & EPH – Phonons, electron–phonon coupling and superconductivity within ABINIT.
7. Molecular & Ab Initio Dynamics
CP2K – BOMD, metadynamics, path integral MD via CP2K/i‑PI interface.
CPMD – Car–Parrinello MD code (plane‑wave pseudopotential), specialized for AIMD.
ABINIT MD – Tutorial and capabilities for BOMD in ABINIT.
i-PI – Python force engine for path integral MD; interfaces to many ab initio codes as clients.
Quantum ESPRESSO CP / PHonon – CPMD and lattice dynamics modules (though CPMD is itself a separate package).
VASP MD – AIMD via VASP (well‑known but not always separately documented; widely used).
8. Structure Prediction & Global Optimization
USPEX – Universal Structure Predictor: Evolutionary Xtallography for crystal structure prediction at variable P–T.
CALYPSO – Crystal structure AnaLYsis by Particle Swarm Optimization (PSO).
AIRSS – Ab initio random structure searching; simple and highly parallel.
XtalOpt – Evolutionary multi‑objective crystal structure prediction (fixed/variable composition).
GASP – Genetic Algorithm for Structure and Phase Prediction; interfaced to VASP, LAMMPS, GULP, etc.
SPINNER – Structure prediction of inorganic crystals using neural network potentials with evolutionary and random search.
EVO – Evolutionary algorithm for crystal structure search.
Minima hopping method – Global optimization via minima hopping; various implementations exist (e.g., in ABINIT/external tools).
MAGUS – Machine learning and graph theory assisted global optimization (mentioned in CSP literature).
9. Post‑Processing, Analysis & Visualization
9.1 DOS, bands, structure analysis
VASPKIT – User‑friendly toolkit for VASP I/O and analysis: elastic, electronic, optical properties; band unfolding; high‑throughput capabilities.
sumo – Python toolkit for plotting/analysis of ab initio solid‑state data; supports VASP, CASTEP, Questaal; DOS, bands, etc.
PyProcar – Python library for PROCAR‑based analysis: band structures, Fermi surfaces (2D/3D), band unfolding, spin/orbital projections.
9.2 Transport (electronic & thermal)
BoltzTraP / BoltzTraP2 – Semi‑classical Boltzmann transport coefficients from band structures; interfaced to WIEN2k, ABINIT, SIESTA, VASP, QE.
AMSET – Efficient electron lifetimes and transport (conductivity, Seebeck, electronic thermal conductivity) from first principles; VASP interface.
EPW (transport module) – Electron transport module in EPW; solves BTE with e–ph scattering; electrical conductivity, mobility, etc.
BoltzWann – Thermoelectric and transport properties using Wannier‑interpolated bands and BTE.
TransOpt – Code for electrical transport of semiconductors (momentum matrix and Boltzmann methods); VASP and QE interfaces.
9.3 Bonding / chemical analysis
LOBSTER – Tool for chemical bonding analysis via projection of PAW DFT wavefunctions onto atomic orbitals; COHP, COOP, charges.
9.4 Fermi surfaces & topology
FermiSurfer – Fermi surface viewer; colors surfaces by Fermi velocity or superconducting gap; 3D Fermi surfaces and sections.
Z2Pack – Calculates topological invariants (Z₂, Chern numbers) via Wannier charge centers; first‑principles and TB modes.
WannierTools – TB analysis and topological classification (Wilson loops, Weyl/Dirac points, Berry curvature).
10. Frameworks, Workflow Engines & Databases
ASE – Atomic Simulation Environment (Python); sets up and steers atomistic simulations; calculators for many codes.
pymatgen – Python Materials Genomics; structures, transformations, analysis, and I/O; powers Materials Project.
AiiDA – Automated Infrastructure and Database for Ab initio Design; workflow engine with provenance tracking; interfaces to QE, VASP, CP2K, etc.
FireWorks – Workflow management engine for scientific workflows; supports dynamic workflows and high‑throughput; used by Materials Project.
atomate – Materials science workflows on top of FireWorks and pymatgen.
atomate2 – Modular workflows for materials science; multi‑code, more flexible than original atomate.
custodian – Just‑in‑time job manager (error detection and recovery) for high‑throughput.
quacc – Python platform for high‑throughput, database‑driven computational materials science and QC.
pyiron – Integrated development environment for computational materials science and workflows.
NOMAD Oasis – Web platform for managing and sharing materials data.
AFLOW – High‑throughput ab initio computing framework (C++).
ASR – Atomic Simulation Recipes, based on ASE.
matminer, pymatflow, httk, matador – Materials informatics toolkits for data mining and workflows.
Materials Project infrastructure – Combines pymatgen, custodian, FireWorks, atomate, databases and APIs.
11. Small, Niche, Community & Research‑Grade Tools
11.1 Niche DFT / ground‑state tools
SeqQuest – Gaussian‑basis pseudopotential DFT (QUEST/SeqQuest).
DACAPO – Plane‑wave pseudopotential DFT code (less widely known).
Socorro – Plane‑wave pseudopotential DFT.
AIMPRO – Localized orbital code for defects in solids.
FreeON (formerly MondoSCF) – Free GTO DFT/HF code.
Firefly (PC GAMESS) – Academic derivative of GAMESS (US) with semi‑empirical capabilities.
11.2 Niche TDDFT / MBPT / excited‑state tools
SALMON – Real‑time TDDFT code for strong fields and electron dynamics (UNCERTAIN placement; primarily known from literature and community lists).
NWChem RT‑TDDFT – Real‑time TDDFT module in NWChem.
Spex TDDFT – TDDFT is implemented but documented as unmaintained; focus is GW/BSE.
11.3 Niche DMFT / many‑body / impurity solvers
ComCTQMC – CT‑QMC impurity solver used in ComDMFT.
ComGutz – Open source LDA+G/Gutzwiller code.
CyGutz – Open source LDA+RISB/Gutzwiller code.
PAUXY – Python‑based auxiliary‑field QMC impurity solvers; phaseless and constrained path AFQMC.
CPMC‑Lab – Matlab package for constrained‑path and phaseless AFQMC on Hubbard models.
11.4 Niche QMC / lattice QMC
StochasticSeriesExpansion.jl – Stochastic series expansion QMC in Julia.
Deterministic QMC (dqmc) – Determinant QMC code for simulating antiferromagnetic quantum critical metals.
11.5 Niche topology / analysis tools
Topological Quantum Chemistry (TQC) – High‑throughput search for topological materials; more of a methodology but associated with code pipelines and databases.
11.6 Niche structure prediction / global optimization
RandSpg – Generates random crystal structures with specified space groups; often used with CSP codes.
MAISE – Module for Ab Initio Structure Evolution; used in structure prediction workflows.
11.7 Niche post‑processing / analysis
API_Phonons – Tool to interface phonons from multiple codes (Phonopy, Phono3py, ShengBTE) and compute thermal transport via Kubo’s linear response.
Crystal Toolkit – Framework for building web apps for materials science; used by Materials Project.
Links and Resources
General lists / overviews
Wikipedia “List of quantum chemistry and solid‑state software” – Baseline list with many codes.
Quest (Sandia) “Alternative DFT codes and web sites” – Lists many periodic and molecular DFT codes.
Psi‑k “Software codes” – List of electronic structure codes for condensed matter.
awesome‑materials‑informatics (GitHub) – Curated list of software, platforms, and datasets in materials informatics.
Ground‑state DFT / MBPT
ABINIT – GW documentation.
CP2K – GW tutorials and documentation.
FHI‑aims – GW/BSE and tutorials.
VASP – Hybrids and GW/BSE/RPA.
Quantum ESPRESSO – What can QE do; hybrids; qe‑acfdt.
Spex – MBPT code for GW/BSE based on FLAPW.
West – GPU‑accelerated GW code.
Yambo – Many‑body code for GW/BSE.
BerkeleyGW – GW and GW‑BSE code.
Fiesta – Gaussian‑basis GW and BSE code.
exciting – FP‑LAPW DFT with excited‑state modules (GW, BSE, TDDFT).
Elk – FP‑LAPW code (general).
DMFT / many‑body
ComDMFT v2.0 paper – Overview of DFT+DMFT and GW+EDMFT codes and impurity solvers (TRIQS, ALPS, iQIST, w2dynamics, EDMFTF, DCore, AMULET, Portobello, etc.).
eDMFT overview – DFT+Embedded DMFT Functional with WIEN2k and CTQMC/OCA/NCA.
DCore paper – Integrated DMFT software built on TRIQS/ALPS.
TRIQS/DFTTools documentation – DFT+DMFT interfaces (Wien2k, Elk, VASP, generic).
ALPS project – ALPS libraries and DMFT impurity solvers; tutorials.
ALPS core libraries paper – Updated core libraries for physics simulations.
Quantum Monte Carlo
QMCPACK manual – AFQMC and real‑space QMC methods.
CASINO review – Overview of variational and diffusion Monte Carlo in CASINO.
TurboRVB website – VMC and LRDMC package.
CHAMP pages – QMC suite (VMC, DMC, optimization).
QWalk – QMC program for molecules and solids.
Tight‑binding / model solvers
Wannier90 website and documentation – MLWFs.
WannierTools docs – Topological analysis in TB framework.
PythTB docs – Tight‑binding and topology package.
TBmodels documentation – Tight‑binding models manipulation and symmetrization.
Tight‑Binding Studio – Tool to construct TB models.
DFTB+ site and docs – Density Functional Tight Binding code.
xtb documentation – Extended tight‑binding (GFN‑xTB) package.
Phonons, lattice dynamics, e–ph
Phonopy documentation – Phonon calculations.
Phono3py documentation – Phonon–phonon interaction and lattice thermal conductivity.
ShengBTE website – BTE for phonons and thermal conductivity.
ALAMODE docs – Anharmonic lattice dynamics.
TDEP – Temperature Dependent Effective Potential code.
EPW documentation – Electron–phonon coupling and transport.
Dynamics
CP2K – Ab initio molecular dynamics.
CPMD – Car–Parrinello MD code.
ABINIT MD tutorial – Molecular dynamics with ABINIT.
i‑PI documentation – Path integral MD interface.
Structure prediction / global optimization
USPEX home and manual – Evolutionary crystal structure prediction.
CALYPSO site and manual – Crystal structure AnaLYsis by Particle Swarm Optimization.
AIRSS documentation – Ab initio random structure searching.
XtalOpt – Evolutionary crystal structure prediction (GitHub/docs).
GASP – Genetic algorithm for structure and phase prediction.
SPINNER – Structure prediction using neural network potentials.
Minima hopping paper – Method and implementations.
Post‑processing / analysis / transport / topology
VASPKIT documentation – VASP I/O and analysis toolkit.
sumo docs – Command‑line tools for solid‑state analysis.
PyProcar docs – PROCAR analysis, Fermi surfaces, band unfolding.
BoltzTraP (TU Wien) – Boltzmann transport properties.
AMSET documentation – Electron transport properties from first principles.
EPW documentation – Electron–phonon coupling and transport.
BoltzWann (Wannier90) – Transport coefficients using Wannier functions.
TransOpt – Electrical transport of semiconductors.
LOBSTER – Bonding analysis from PAW DFT.
FermiSurfer – Fermi surface visualization.
Z2Pack – Topological invariants calculation.
Frameworks / workflows / databases
ASE – Atomic Simulation Environment.
pymatgen – Python Materials Genomics.
AiiDA – Automated infrastructure and database for ab initio design.
FireWorks – Workflow engine developed at LBNL.
atomate and atomate2 – Materials science workflows based on pymatgen/FireWorks.
custodian – Job management framework for high‑throughput.
awesome‑materials‑informatics – Curated list including quacc, pyiron, NOMAD Oasis, etc.


mindmap
  root((Computational Tools))
    Ground-State DFT and Variants
      VASP, Quantum ESPRESSO, ABINIT, CASTEP
      CP2K, GPAW, SIESTA, OpenMX, FHI-aims
      WIEN2k, FLEUR, ELK, exciting, RSPt
      CRYSTAL, ORCA, Q-Chem, Gaussian, NWChem, PSI4, PySCF
    Time-Dependent and Excited-State
      Octopus, GPAW RT-TDDFT
      Yambo, BerkeleyGW, VASP GW/BSE
      Spex, exciting, West, FHI-aims GW/BSE
      Fiesta, FHI-aims TDDFT
    Strongly Correlated / Many-Body
      TRIQS, DCore, EDMFTF, ComDMFT, AMULET
      DMFTwDFT, w2dynamics, iQIST, ALPS CT-QMC
      QMCPACK, CASINO, TurboRVB, CHAMP, QWalk
    Wavefunction-Based QC
      ORCA, CFOUR, MRCC, PSI4, PySCF
      OpenMolcas, Molpro, Dalton, DIRAC, COLUMBUS, eT
    Tight-Binding and Model Solvers
      Wannier90, WannierTools, PythTB, TBmodels
      TRIQS/DFTTools
    Phonons, Lattice Dynamics, e–ph
      Phonopy, phono3py, ShengBTE, almaBTE
      ALAMODE, TDEP, EPW, ABINIT DFPT, QE phonons
    Molecular and Ab Initio Dynamics
      CP2K AIMD, CPMD, ABINIT MD, i-PI
    Structure Prediction / Global Optimization
      USPEX, CALYPSO, AIRSS, XtalOpt, GASP
      SPINNER, EVO, Minima Hopping (various impls)
    Post-Processing / Analysis
      VASPKIT, LOBSTER, BoltzTraP, BoltzTraP2
      EPW transport, AMSET, BoltzWann, TransOpt
      sumo, PyProcar, FermiSurfer, Z2Pack
    Frameworks / Workflows / Databases
      ASE, pymatgen, AiiDA, FireWorks, atomate, atomate2
      custodian, quacc, pyiron, NOMAD Oasis, AFLOW
    Small / Niche / Research-Grade
      Many specialized TB, DMFT, QMC, topology codes (listed below)