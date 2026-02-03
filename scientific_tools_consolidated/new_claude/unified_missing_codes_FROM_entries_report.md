# Unified vs entries*.md: Missing Codes Report

- Unified file: `scientific_tools_consolidated/new_claude/unified_all_codes_VERIFIED_ENHANCED.md`
- Entries files glob: `scientific_tools_consolidated/new_claude/entries*.md`

## Summary
- Unique code names in entries*.md: **279**
- Unique code names in unified file: **279**
- Codes present in entries*.md but **missing in unified**: **64**
- Codes present in unified but **not present in entries*.md**: **64**

## Missing-vs-Rename Heuristic
- Missing codes with obvious name-overlap candidate(s) in unified: **23 / 64**

## Missing In Unified (present in entries*.md)

### ABINIT-GW
- IDs in entries tables: 50
- Context (nearest headings in source): Category: 1.1 PLANE-WAVE PSEUDOPOTENTIAL & PAW METHODS > Category: 1.2 ALL-ELECTRON & FULL-POTENTIAL METHODS > Category: 1.3 QUANTUM CHEMISTRY (GAUSSIAN BASIS)
- Source rows: `scientific_tools_consolidated/new_claude/entries_001_100_complete.md:70`
- License: GNU GPL v2
- Official website: Part of ABINIT
- Basis / method: GW
- Specialization: GW, BSE in ABINIT
- Primary publication (as in entries table): Gonze, X.; et al. "ABINIT: First-principles approach." *Comput. Phys. Commun.* **180**, 2582-2615 (2009). https://doi.org/10.1016/j.cpc.2009.02.005
- Unified name-overlap candidate(s): ABINIT
- Unified closest-by-similarity (normalized): ABINIT, Pybinding

### ACE
- IDs in entries tables: 320
- Context (nearest headings in source): Category: 11.6 SPECIALIZED SOLVERS & METHODS > Category: 11.7 ADDITIONAL SPECIALIZED CODES > Category: 11.7.2 ADDITIONAL FRAMEWORKS & UTILITIES
- Source rows: `scientific_tools_consolidated/new_claude/entries_291_320_complete.md:47`
- License: MIT License
- Official website: https://github.com/ACEsuit/ACE.jl
- Basis / method: Equivariant ML potentials
- Specialization: Atomic Cluster Expansion, complete basis, body-ordered descriptors, Julia implementation
- Primary publication (as in entries table): Drautz, R. "Atomic cluster expansion for accurate and transferable interatomic potentials." *Phys. Rev. B* **99**, 014104 (2019). https://doi.org/10.1103/PhysRevB.99.014104; Dusson, G.; Bachmayr, M.; Csányi, G.; Drautz, R.; Etter, S.; van der Oord, C.; Ortner, C. "Atomic cluster expansion: Completeness, efficiency and stability." *J. Comput. Phys.* **454**, 110946 (2022). https://doi.org/10.1016/j.jcp.2022.110946
- Unified name-overlap candidate(s): ACES II, MACE
- Unified closest-by-similarity (normalized): MACE, ASE, ADC, ACES II

### ACES
- IDs in entries tables: 92
- Context (nearest headings in source): Category: 1.1 PLANE-WAVE PSEUDOPOTENTIAL & PAW METHODS > Category: 1.2 ALL-ELECTRON & FULL-POTENTIAL METHODS > Category: 1.3 QUANTUM CHEMISTRY (GAUSSIAN BASIS)
- Source rows: `scientific_tools_consolidated/new_claude/entries_001_100_complete.md:121`
- License: Academic
- Official website: http://www.qtp.ufl.edu/Aces/
- Basis / method: CC methods
- Specialization: Advanced Coupled Electron Systems
- Primary publication (as in entries table): Stanton, J. F. "ACES II coupled-cluster package." (1991+). http://www.qtp.ufl.edu/Aces/
- Unified name-overlap candidate(s): ACES II
- Unified closest-by-similarity (normalized): ACES II, MACE, ASE, ADC

### ADF-BAND
- IDs in entries tables: 114
- Context (nearest headings in source): Category: 1.5 ALL-ELECTRON NUMERIC BASIS > Category: 1.6 LINEAR-SCALING & ORDER-N METHODS > Category: 1.7 SPECIALIZED DFT IMPLEMENTATIONS
- Source rows: `scientific_tools_consolidated/new_claude/entries_101_200_complete.md:36`
- License: Commercial
- Official website: https://www.scm.com/product/band/
- Basis / method: STOs (periodic)
- Specialization: Periodic STO basis, ADF suite
- Primary publication (as in entries table): te Velde, G.; Baerends, E. J. "Precise density-functional method for periodic structures." *Phys. Rev. B* **44**, 7888-7903 (1991). https://doi.org/10.1103/PhysRevB.44.7888
- Unified name-overlap candidate(s): ADF
- Unified closest-by-similarity (normalized): BandUP, ADF

### Amber
- IDs in entries tables: 171
- Context (nearest headings in source): Category: 4. TIGHT-BINDING & WANNIER FUNCTIONS > Category: 5. PHONONS & LATTICE DYNAMICS > Category: 6. MOLECULAR DYNAMICS
- Source rows: `scientific_tools_consolidated/new_claude/entries_101_200_complete.md:123`
- License: Commercial/Academic
- Official website: https://ambermd.org/
- Basis / method: Classical MD (biomolecular)
- Specialization: Biomolecular MD, force fields
- Primary publication (as in entries table): Case, D. A.; et al. "The Amber biomolecular simulation programs." *J. Comput. Chem.* **26**, 1668-1688 (2005). https://doi.org/10.1002/jcc.20290
- Unified closest-by-similarity (normalized): almaBTE, matminer, Yambo, FLAME, Bader

### Atompaw
- IDs in entries tables: 299, 336
- Context (nearest headings in source): COMPREHENSIVE MATERIALS CODES CATALOG - ENTRIES 321-400 > Triple-Verified with 100% Accuracy Standard | January 2026 > Category: 9.1 ELECTRONIC ANALYSIS & POST-PROCESSING (Continued); Triple-Verified with 100% Accuracy Standard | January 2026 > Category: 11.6 SPECIALIZED SOLVERS & METHODS > Category: 11.7 ADDITIONAL SPECIALIZED CODES
- Source rows: `scientific_tools_consolidated/new_claude/entries_291_320_complete.md:21`, `scientific_tools_consolidated/new_claude/entries_321_400_complete.md:23`
- License: GNU GPL
- Official website: https://github.com/AtomPAW/atompaw
- Basis / method: PAW dataset generation; PAW generator
- Specialization: PAW pseudopotential generator, ABINIT/PWSCF compatible, all-electron reference; PAW pseudopotential generator, all-electron
- Primary publication (as in entries table): Holzwarth, N. A. W.; Tackett, A. R.; Matthews, G. E. "A Projector Augmented Wave (PAW) code for electronic structure calculations, Part I: atompaw for generating atom-centered functions." *Comput. Phys. Commun.* **135**, 329-347 (2001). https://doi.org/10.1016/S0010-4655(00)00244-7; Holzwarth, N. A. W.; et al. "A Projector Augmented Wave (PAW) code." *Comput. Phys. Commun.* **135**, 329-347 (2001). https://doi.org/10.1016/S0010-4655(00)00244-7
- Unified name-overlap candidate(s): Atompaw/PWPAW
- Unified closest-by-similarity (normalized): Atompaw/PWPAW, atomate, atomate2, FSatom, CatMAP

### AtomPAW
- IDs in entries tables: 84
- Context (nearest headings in source): Category: 1.1 PLANE-WAVE PSEUDOPOTENTIAL & PAW METHODS > Category: 1.2 ALL-ELECTRON & FULL-POTENTIAL METHODS > Category: 1.3 QUANTUM CHEMISTRY (GAUSSIAN BASIS)
- Source rows: `scientific_tools_consolidated/new_claude/entries_001_100_complete.md:113`
- License: GNU GPL
- Official website: https://github.com/AtomPAW/atompaw
- Basis / method: PAW generator
- Specialization: PAW pseudopotential generator
- Primary publication (as in entries table): Holzwarth, N. A. W.; et al. "A Projector Augmented Wave (PAW) code." *Comput. Phys. Commun.* **135**, 329-347 (2001). https://doi.org/10.1016/S0010-4655(00)00244-7
- Unified name-overlap candidate(s): Atompaw/PWPAW
- Unified closest-by-similarity (normalized): Atompaw/PWPAW, atomate, atomate2, FSatom, CatMAP

### Block
- IDs in entries tables: 134
- Context (nearest headings in source): Category: 2. EXCITED-STATE METHODS (GW/BSE/TDDFT) > Category: 3. STRONGLY CORRELATED SYSTEMS > Category: 3.2 MODEL HAMILTONIANS & PEDAGOGICAL
- Source rows: `scientific_tools_consolidated/new_claude/entries_101_200_complete.md:71`
- License: Research
- Official website: https://sanshar.github.io/Block/
- Basis / method: DMRG
- Specialization: DMRG for quantum chemistry
- Primary publication (as in entries table): Sharma, S.; Chan, G. K.-L. "Spin-adapted density matrix renormalization group algorithms." *J. Chem. Phys.* **136**, 124121 (2012). https://doi.org/10.1063/1.3695642

### BOSS
- IDs in entries tables: 186
- Context (nearest headings in source): Category: 6. MOLECULAR DYNAMICS > Category: 7. TRANSITION STATES & RARE EVENTS > Category: 8. STRUCTURE PREDICTION
- Source rows: `scientific_tools_consolidated/new_claude/entries_101_200_complete.md:148`
- License: Research
- Official website: Research code
- Basis / method: Basin hopping
- Specialization: Basin hopping, global optimization
- Primary publication (as in entries table): Wales, D. J.; Doye, J. P. K. "Global Optimization by Basin-Hopping." *J. Phys. Chem. A* **101**, 5111-5116 (1997). https://doi.org/10.1021/jp970984n

### BSE
- IDs in entries tables: 88
- Context (nearest headings in source): Category: 1.1 PLANE-WAVE PSEUDOPOTENTIAL & PAW METHODS > Category: 1.2 ALL-ELECTRON & FULL-POTENTIAL METHODS > Category: 1.3 QUANTUM CHEMISTRY (GAUSSIAN BASIS)
- Source rows: `scientific_tools_consolidated/new_claude/entries_001_100_complete.md:117`
- License: Academic
- Official website: https://www.basissetexchange.org/
- Basis / method: Basis sets
- Specialization: Basis Set Exchange, repository
- Primary publication (as in entries table): Pritchard, B. P.; et al. "New basis set exchange." *J. Chem. Inf. Model.* **59**, 4814-4820 (2019). https://doi.org/10.1021/acs.jcim.9b00725
- Unified name-overlap candidate(s): NBSE
- Unified closest-by-similarity (normalized): NBSE, ASE, Lobster, Spex

### CCDC
- IDs in entries tables: 89
- Context (nearest headings in source): Category: 1.1 PLANE-WAVE PSEUDOPOTENTIAL & PAW METHODS > Category: 1.2 ALL-ELECTRON & FULL-POTENTIAL METHODS > Category: 1.3 QUANTUM CHEMISTRY (GAUSSIAN BASIS)
- Source rows: `scientific_tools_consolidated/new_claude/entries_001_100_complete.md:118`
- License: Commercial/Academic
- Official website: https://www.ccdc.cam.ac.uk/
- Basis / method: Crystal database
- Specialization: Cambridge Structural Database, experimental
- Primary publication (as in entries table): Groom, C. R.; et al. "The Cambridge Structural Database." *Acta Crystallogr. B* **72**, 171-179 (2016). https://doi.org/10.1107/S2052520616003954
- Unified closest-by-similarity (normalized): DCA++, ADC

### CHARMM
- IDs in entries tables: 172
- Context (nearest headings in source): Category: 4. TIGHT-BINDING & WANNIER FUNCTIONS > Category: 5. PHONONS & LATTICE DYNAMICS > Category: 6. MOLECULAR DYNAMICS
- Source rows: `scientific_tools_consolidated/new_claude/entries_101_200_complete.md:124`
- License: Academic/Commercial
- Official website: https://www.charmm.org/
- Basis / method: Classical MD
- Specialization: Chemistry at Harvard Macromolecular Mechanics
- Primary publication (as in entries table): Brooks, B. R.; et al. "CHARMM: The Biomolecular Simulation Program." *J. Comput. Chem.* **30**, 1545-1614 (2009). https://doi.org/10.1002/jcc.21287
- Unified closest-by-similarity (normalized): CHAMP

### CHGNet
- IDs in entries tables: 280
- Context (nearest headings in source): Category: 10.2 HIGH-THROUGHPUT & DATABASE INFRASTRUCTURE > Category: 10.2.2 SPECIALIZED HIGH-THROUGHPUT TOOLS > Category: 11.3 MACHINE LEARNING POTENTIALS & NEURAL NETWORKS
- Source rows: `scientific_tools_consolidated/new_claude/entries_261_290_complete.md:37`
- License: BSD 3-Clause
- Official website: https://github.com/CederGroupHub/chgnet
- Basis / method: Crystal Hamiltonian GNN
- Specialization: Charge-informed, magnetic moments, forces & stresses
- Primary publication (as in entries table): Deng, B.; Zhong, P.; Jun, K.; Riebesell, J.; Han, K.; Bartel, C. J.; Ceder, G. "CHGNet as a pretrained universal neural network potential for charge-informed atomistic modelling." *Nat. Mach. Intell.* **5**, 1031-1041 (2023). https://doi.org/10.1038/s42256-023-00716-3
- Unified closest-by-similarity (normalized): SchNetPack, CONQUEST

### COPEX
- IDs in entries tables: 189
- Context (nearest headings in source): Category: 6. MOLECULAR DYNAMICS > Category: 7. TRANSITION STATES & RARE EVENTS > Category: 8. STRUCTURE PREDICTION
- Source rows: `scientific_tools_consolidated/new_claude/entries_101_200_complete.md:151`
- License: Research
- Official website: Research code
- Basis / method: Co-evolutionary
- Specialization: Co-evolutionary for complex systems
- Primary publication (as in entries table): Wang, Y.; et al. "COPEX: co-evolutionary crystal structure prediction algorithm." *npj Comput. Mater.* **7**, 195 (2021). https://doi.org/10.1038/s41524-021-00668-5
- Unified closest-by-similarity (normalized): OpenMX, Spex, COHP, USPEX, DCore

### CPMD
- IDs in entries tables: 168, 233
- Context (nearest headings in source): Category: 4. TIGHT-BINDING & WANNIER FUNCTIONS > Category: 5. PHONONS & LATTICE DYNAMICS > Category: 6. MOLECULAR DYNAMICS; Section 7: Molecular & Ab Initio Dynamics > Verified with 100% Accuracy Standard | January 2026 > Category: 7.1 Born-Oppenheimer Molecular Dynamics & Path Integral Methods
- Source rows: `scientific_tools_consolidated/new_claude/entries_101_200_complete.md:120`, `scientific_tools_consolidated/new_claude/entries_201_260_complete.md:169`, `scientific_tools_consolidated/new_claude/entries_231_233_molecular_dynamics.md:11`
- License: MIT License (2022+)
- Official website: https://github.com/CPMD-code<br>Historical: http://www.cpmd.org/; https://www.cpmd.org/
- Basis / method: Car-Parrinello MD; Plane-wave/pseudopotential AIMD
- Specialization: Car-Parrinello MD (fictitious electronic dynamics), plane-wave DFT, QM/MM with GROMOS, path integral MD, metadynamics; Car-Parrinello molecular dynamics
- Primary publication (as in entries table): Car, R.; Parrinello, M. "Unified Approach for Molecular Dynamics and Density-Functional Theory." *Phys. Rev. Lett.* **55** (22), 2471-2474 (1985). https://doi.org/10.1103/PhysRevLett.55.2471<br><br>Hutter, J.; Iannuzzi, M.; Schiffmann, F.; VandeVondele, J. "cp2k: atomistic simulations of condensed matter systems." *Wiley Interdiscip. Rev.: Comput. Mol. Sci.* **4** (1), 15-25 (2014). https://doi.org/10.1002/wcms.1159; Car, R.; Parrinello, M. "Unified Approach for Molecular Dynamics and Density-Functional Theory." *Phys. Rev. Lett.* **55**, 2471-2474 (1985). https://doi.org/10.1103/PhysRevLett.55.2471
- Unified closest-by-similarity (normalized): PLUMED, VMD

### DL_POLY
- IDs in entries tables: 174
- Context (nearest headings in source): Category: 4. TIGHT-BINDING & WANNIER FUNCTIONS > Category: 5. PHONONS & LATTICE DYNAMICS > Category: 6. MOLECULAR DYNAMICS
- Source rows: `scientific_tools_consolidated/new_claude/entries_101_200_complete.md:126`
- License: BSD 3-Clause
- Official website: https://www.scd.stfc.ac.uk/Pages/DL_POLY.aspx
- Basis / method: Classical MD
- Specialization: Daresbury Laboratory polymer MD
- Primary publication (as in entries table): Todorov, I. T.; et al. "DL_POLY_3: new dimensions in molecular dynamics simulations." *J. Mater. Chem.* **16**, 1911-1918 (2006). https://doi.org/10.1039/B517931A

### DMDW
- IDs in entries tables: 165, 230
- Context (nearest headings in source): Category: 3.2 MODEL HAMILTONIANS & PEDAGOGICAL > Category: 4. TIGHT-BINDING & WANNIER FUNCTIONS > Category: 5. PHONONS & LATTICE DYNAMICS; Category: 4. WAVEFUNCTION-BASED QUANTUM CHEMISTRY > Category: 5. TIGHT-BINDING, MODEL HAMILTONIANS & DOWNFOLDING > Category: 6. PHONONS, LATTICE DYNAMICS & ELECTRON-PHONON
- Source rows: `scientific_tools_consolidated/new_claude/entries_101_200_complete.md:112`, `scientific_tools_consolidated/new_claude/entries_201_230_materials_codes.md:47`, `scientific_tools_consolidated/new_claude/entries_201_260_complete.md:47`
- License: Academic; Research
- Official website: Research code; http://www.dmdw.de/
- Basis / method: Debye-Waller factors; Dynamic matrix
- Specialization: Debye-Waller, thermal motion, X-ray; Dynamic matrix approaches
- Primary publication (as in entries table): Dynamic matrix and phonon dispersion calculations; Kresse, G.; Furthmüller, J.; Hafner, J. "Ab initio Force Constant Approach to Phonon Dispersion Relations of Diamond and Graphite." *Europhys. Lett.* **32**, 729 (1995). https://doi.org/10.1209/0295-5075/32/9/005 (DMDW methodology for phonon properties)
- Unified name-overlap candidate(s): DMDW/RTDW
- Unified closest-by-similarity (normalized): DMDW/RTDW, VMD

### DMOL3
- IDs in entries tables: 106
- Context (nearest headings in source): Triple-Verified with 100% Accuracy Standard | January 2026 > Category: 1.4 GAUSSIAN BASIS - SPECIALIZED PACKAGES > Category: 1.5 ALL-ELECTRON NUMERIC BASIS
- Source rows: `scientific_tools_consolidated/new_claude/entries_101_200_complete.md:18`
- License: Commercial
- Official website: Part of BIOVIA
- Basis / method: Numeric basis
- Specialization: Numeric atomic orbitals, BIOVIA Materials Studio
- Primary publication (as in entries table): Delley, B. "From molecules to solids with the DMol3 approach." *J. Chem. Phys.* **113**, 7756-7764 (2000). https://doi.org/10.1063/1.1316015
- Unified closest-by-similarity (normalized): JMol, molgw, PyMOL

### DMRG++
- IDs in entries tables: 133
- Context (nearest headings in source): Category: 2. EXCITED-STATE METHODS (GW/BSE/TDDFT) > Category: 3. STRONGLY CORRELATED SYSTEMS > Category: 3.2 MODEL HAMILTONIANS & PEDAGOGICAL
- Source rows: `scientific_tools_consolidated/new_claude/entries_101_200_complete.md:70`
- License: BSD 2-Clause
- Official website: https://g1257.github.io/dmrgPlusPlus/
- Basis / method: DMRG
- Specialization: Density Matrix Renormalization Group, Oak Ridge
- Primary publication (as in entries table): Alvarez, G. "The density matrix renormalization group for strongly correlated electron systems." *Comput. Phys. Commun.* **180**, 1572-1578 (2009). https://doi.org/10.1016/j.cpc.2009.02.016
- Unified closest-by-similarity (normalized): RMG

### Dual fermions
- IDs in entries tables: 353
- Context (nearest headings in source): Category: 9.1 ELECTRONIC ANALYSIS & POST-PROCESSING (Continued) > Category: 9.2 TRANSPORT & OPTICAL PROPERTIES > Category: 9.3 SPECIALIZED DFT CODES (Additional)
- Source rows: `scientific_tools_consolidated/new_claude/entries_321_400_complete.md:50`
- License: Research
- Official website: Various implementations
- Basis / method: Beyond DMFT
- Specialization: Diagrammatic extension of DMFT
- Primary publication (as in entries table): Rubtsov, A. N.; et al. "Dual fermion approach to nonlocal correlations." *Phys. Rev. B* **77**, 033101 (2008). https://doi.org/10.1103/PhysRevB.77.033101
- Unified name-overlap candidate(s): dual fermions
- Unified closest-by-similarity (normalized): dual fermions, SALMON, Dalton

### EDLib
- IDs in entries tables: 139
- Context (nearest headings in source): Category: 2. EXCITED-STATE METHODS (GW/BSE/TDDFT) > Category: 3. STRONGLY CORRELATED SYSTEMS > Category: 3.2 MODEL HAMILTONIANS & PEDAGOGICAL
- Source rows: `scientific_tools_consolidated/new_claude/entries_101_200_complete.md:76`
- License: MIT License
- Official website: https://github.com/Q-solvers/EDLib
- Basis / method: Exact diagonalization
- Specialization: Exact diagonalization for impurity models
- Primary publication (as in entries table): Iskakov, S.; et al. "EDLib: Exact diagonalization library." (GitHub). https://github.com/Q-solvers/EDLib

### EXC
- IDs in entries tables: 125
- Context (nearest headings in source): Category: 1.6 LINEAR-SCALING & ORDER-N METHODS > Category: 1.7 SPECIALIZED DFT IMPLEMENTATIONS > Category: 2. EXCITED-STATE METHODS (GW/BSE/TDDFT)
- Source rows: `scientific_tools_consolidated/new_claude/entries_101_200_complete.md:52`
- License: Part of exciting
- Official website: Part of exciting
- Basis / method: BSE/TDDFT
- Specialization: Excited states in exciting
- Primary publication (as in entries table): Gulans, A.; et al. (See exciting documentation)
- Unified name-overlap candidate(s): exciting
- Unified closest-by-similarity (normalized): EMC, Spex, DDEC

### FHI-gap
- IDs in entries tables: 48
- Context (nearest headings in source): Category: 1.1 PLANE-WAVE PSEUDOPOTENTIAL & PAW METHODS > Category: 1.2 ALL-ELECTRON & FULL-POTENTIAL METHODS > Category: 1.3 QUANTUM CHEMISTRY (GAUSSIAN BASIS)
- Source rows: `scientific_tools_consolidated/new_claude/entries_001_100_complete.md:68`
- License: Academic
- Official website: Part of FHI-aims
- Basis / method: GW/RPA
- Specialization: GW, RPA, MP2 in numeric AO basis
- Primary publication (as in entries table): Ren, X.; et al. "Resolution-of-identity approach to Hartree–Fock, hybrid density functionals, RPA, MP2 and GW." *New J. Phys.* **14**, 053020 (2012). https://doi.org/10.1088/1367-2630/14/5/053020
- Unified closest-by-similarity (normalized): HiLAPW, FHI-aims, GASP

### Fleur-inpgen
- IDs in entries tables: 97
- Context (nearest headings in source): Category: 1.1 PLANE-WAVE PSEUDOPOTENTIAL & PAW METHODS > Category: 1.2 ALL-ELECTRON & FULL-POTENTIAL METHODS > Category: 1.3 QUANTUM CHEMISTRY (GAUSSIAN BASIS)
- Source rows: `scientific_tools_consolidated/new_claude/entries_001_100_complete.md:126`
- License: Academic
- Official website: https://www.flapw.de/
- Basis / method: LAPW utilities
- Specialization: Fleur input generator
- Primary publication (as in entries table): Part of Fleur package for input generation
- Unified name-overlap candidate(s): Fleur
- Unified closest-by-similarity (normalized): Fleur

### GOLLUM
- IDs in entries tables: 119
- Context (nearest headings in source): Category: 1.5 ALL-ELECTRON NUMERIC BASIS > Category: 1.6 LINEAR-SCALING & ORDER-N METHODS > Category: 1.7 SPECIALIZED DFT IMPLEMENTATIONS
- Source rows: `scientific_tools_consolidated/new_claude/entries_101_200_complete.md:41`
- License: Academic
- Official website: https://www.lancaster.ac.uk/physics/research/condensed-matter/quantum-technology/gollum/
- Basis / method: Quantum transport
- Specialization: Transport, thermoelectric, spin
- Primary publication (as in entries table): Ferrer, J.; et al. "GOLLUM: a next-generation simulation tool for electron, thermal and spin transport." *New J. Phys.* **16**, 093029 (2014). https://doi.org/10.1088/1367-2630/16/9/093029
- Unified closest-by-similarity (normalized): Columbus

### GROMACS
- IDs in entries tables: 169
- Context (nearest headings in source): Category: 4. TIGHT-BINDING & WANNIER FUNCTIONS > Category: 5. PHONONS & LATTICE DYNAMICS > Category: 6. MOLECULAR DYNAMICS
- Source rows: `scientific_tools_consolidated/new_claude/entries_101_200_complete.md:121`
- License: GNU LGPL v2.1
- Official website: https://www.gromacs.org/
- Basis / method: Classical MD
- Specialization: Biomolecular MD specialist
- Primary publication (as in entries table): Abraham, M. J.; et al. "GROMACS: High performance molecular simulations." *SoftwareX* **1-2**, 19-25 (2015). https://doi.org/10.1016/j.softx.2015.06.001

### GW4
- IDs in entries tables: 46
- Context (nearest headings in source): Category: 1.1 PLANE-WAVE PSEUDOPOTENTIAL & PAW METHODS > Category: 1.2 ALL-ELECTRON & FULL-POTENTIAL METHODS > Category: 1.3 QUANTUM CHEMISTRY (GAUSSIAN BASIS)
- Source rows: `scientific_tools_consolidated/new_claude/entries_001_100_complete.md:66`
- License: Academic
- Official website: Research code
- Basis / method: GW
- Specialization: GW self-energy calculations
- Primary publication (as in entries table): Various implementations of GW approximation
- Unified closest-by-similarity (normalized): GPAW

### HOOMD-blue
- IDs in entries tables: 173
- Context (nearest headings in source): Category: 4. TIGHT-BINDING & WANNIER FUNCTIONS > Category: 5. PHONONS & LATTICE DYNAMICS > Category: 6. MOLECULAR DYNAMICS
- Source rows: `scientific_tools_consolidated/new_claude/entries_101_200_complete.md:125`
- License: BSD 3-Clause
- Official website: https://glotzerlab.engin.umich.edu/hoomd-blue/
- Basis / method: Particle MD
- Specialization: GPU-accelerated particle simulations
- Primary publication (as in entries table): Anderson, J. A.; et al. "HOOMD-blue: A Python package for high-performance molecular dynamics." *Comput. Mater. Sci.* **173**, 109363 (2020). https://doi.org/10.1016/j.commatsci.2019.109363

### HORTON
- IDs in entries tables: 99
- Context (nearest headings in source): Category: 1.1 PLANE-WAVE PSEUDOPOTENTIAL & PAW METHODS > Category: 1.2 ALL-ELECTRON & FULL-POTENTIAL METHODS > Category: 1.3 QUANTUM CHEMISTRY (GAUSSIAN BASIS)
- Source rows: `scientific_tools_consolidated/new_claude/entries_001_100_complete.md:128`
- License: GNU GPL v3
- Official website: http://theochem.github.io/horton/
- Basis / method: QC development
- Specialization: Modular QC development platform
- Primary publication (as in entries table): Verstraelen, T.; et al. "HORTON: Helpful Open-source Research TOol." (2011-2018). http://theochem.github.io/horton/
- Unified closest-by-similarity (normalized): PHONON, PHON

### HTST
- IDs in entries tables: 180
- Context (nearest headings in source): Category: 5. PHONONS & LATTICE DYNAMICS > Category: 6. MOLECULAR DYNAMICS > Category: 7. TRANSITION STATES & RARE EVENTS
- Source rows: `scientific_tools_consolidated/new_claude/entries_101_200_complete.md:137`
- License: Research
- Official website: Research codes
- Basis / method: Harmonic TST
- Specialization: Harmonic Transition State Theory
- Primary publication (as in entries table): Henkelman, G.; Jónsson, H. "Long time scale kinetic Monte Carlo simulations." *J. Chem. Phys.* **115**, 9657-9666 (2001). https://doi.org/10.1063/1.1415500
- Unified closest-by-similarity (normalized): PhonTS, HTOCSP, HOTBIT

### Hubbard-I
- IDs in entries tables: 129
- Context (nearest headings in source): Category: 1.7 SPECIALIZED DFT IMPLEMENTATIONS > Category: 2. EXCITED-STATE METHODS (GW/BSE/TDDFT) > Category: 3. STRONGLY CORRELATED SYSTEMS
- Source rows: `scientific_tools_consolidated/new_claude/entries_101_200_complete.md:61`
- License: Part of TRIQS
- Official website: Hubbard-I approximation
- Basis / method: Hubbard, J. "Electron correlations in narrow energy bands." *Proc. R. Soc. London Ser. A* **276**, 238-257 (1963). https://doi.org/10.1098/rspa.1963.0204
- Primary publication (as in entries table): Atomic limit DMFT solver
- Unified closest-by-similarity (normalized): HubbardFermiMatsubara

### iTensor
- IDs in entries tables: 131
- Context (nearest headings in source): Category: 2. EXCITED-STATE METHODS (GW/BSE/TDDFT) > Category: 3. STRONGLY CORRELATED SYSTEMS > Category: 3.2 MODEL HAMILTONIANS & PEDAGOGICAL
- Source rows: `scientific_tools_consolidated/new_claude/entries_101_200_complete.md:68`
- License: Apache 2.0
- Official website: https://itensor.org/
- Basis / method: Tensor networks
- Specialization: DMRG, MPS, tensor networks, C++/Julia
- Primary publication (as in entries table): Fishman, M.; et al. "The ITensor Software Library for Tensor Network Calculations." *SciPost Phys. Codebases* **4** (2022). https://doi.org/10.21468/SciPostPhysCodeb.4

### Julia packages
- IDs in entries tables: 85
- Context (nearest headings in source): Category: 1.1 PLANE-WAVE PSEUDOPOTENTIAL & PAW METHODS > Category: 1.2 ALL-ELECTRON & FULL-POTENTIAL METHODS > Category: 1.3 QUANTUM CHEMISTRY (GAUSSIAN BASIS)
- Source rows: `scientific_tools_consolidated/new_claude/entries_001_100_complete.md:114`
- License: MIT
- Official website: https://github.com/JuliaMolSim
- Basis / method: Julia ecosystem
- Specialization: Julia-based materials codes
- Primary publication (as in entries table): Stukowski, A. "Computational Materials Science with Julia." Various GitHub repos. https://github.com/JuliaMolSim
- Unified closest-by-similarity (normalized): Julia materials packages

### Libint
- IDs in entries tables: 100
- Context (nearest headings in source): Category: 1.1 PLANE-WAVE PSEUDOPOTENTIAL & PAW METHODS > Category: 1.2 ALL-ELECTRON & FULL-POTENTIAL METHODS > Category: 1.3 QUANTUM CHEMISTRY (GAUSSIAN BASIS)
- Source rows: `scientific_tools_consolidated/new_claude/entries_001_100_complete.md:129`
- License: GNU LGPL v3
- Official website: https://github.com/evaleev/libint
- Basis / method: Integral library
- Specialization: Gaussian integral evaluation library
- Primary publication (as in entries table): Valeev, E. F. "Libint: A library for the evaluation of molecular integrals." https://github.com/evaleev/libint
- Unified closest-by-similarity (normalized): ABINIT

### Magnon codes
- IDs in entries tables: 378
- Context (nearest headings in source): Category: 9.3 SPECIALIZED DFT CODES (Additional) > Category: 10.1 WORKFLOW ENGINES & FRAMEWORKS > Category: 10.2 DATABASES & HIGH-THROUGHPUT
- Source rows: `scientific_tools_consolidated/new_claude/entries_321_400_complete.md:85`
- License: Various
- Official website: Various
- Basis / method: Magnon dynamics
- Specialization: Spin wave/magnon calculations
- Primary publication (as in entries table): Various implementations for magnonic properties
- Unified name-overlap candidate(s): magnon codes
- Unified closest-by-similarity (normalized): magnon codes

### Matgl
- IDs in entries tables: 279
- Context (nearest headings in source): Category: 10.2 HIGH-THROUGHPUT & DATABASE INFRASTRUCTURE > Category: 10.2.2 SPECIALIZED HIGH-THROUGHPUT TOOLS > Category: 11.3 MACHINE LEARNING POTENTIALS & NEURAL NETWORKS
- Source rows: `scientific_tools_consolidated/new_claude/entries_261_290_complete.md:36`
- License: BSD 3-Clause
- Official website: https://github.com/materialsvirtuallab/matgl
- Basis / method: Graph neural networks
- Specialization: M3GNet, megaNet, universal periodic table coverage
- Primary publication (as in entries table): Chen, C.; Ong, S. P. "A universal graph deep learning interatomic potential for the periodic table." *Nat. Comput. Sci.* **2**, 718-728 (2022). https://doi.org/10.1038/s43588-022-00349-3
- Unified closest-by-similarity (normalized): pymatgen, MatPy, BAGEL

### MLIP
- IDs in entries tables: 276, 374
- Context (nearest headings in source): Category: 10.2 HIGH-THROUGHPUT & DATABASE INFRASTRUCTURE > Category: 10.2.2 SPECIALIZED HIGH-THROUGHPUT TOOLS > Category: 11.3 MACHINE LEARNING POTENTIALS & NEURAL NETWORKS; Category: 9.3 SPECIALIZED DFT CODES (Additional) > Category: 10.1 WORKFLOW ENGINES & FRAMEWORKS > Category: 10.2 DATABASES & HIGH-THROUGHPUT
- Source rows: `scientific_tools_consolidated/new_claude/entries_261_290_complete.md:33`, `scientific_tools_consolidated/new_claude/entries_321_400_complete.md:81`
- License: MIT License
- Official website: https://gitlab.com/ashapeev/mlip-2; https://gitlab.com/ashapeev/mlip-3
- Basis / method: ML potentials; Moment tensor potentials
- Specialization: MTP method, active learning, linear complexity; Moment tensor potentials, Skoltech
- Primary publication (as in entries table): Shapeev, A. V. "Moment Tensor Potentials: A class of systematically improvable interatomic potentials." *Multiscale Model. Simul.* **14**, 1153-1173 (2016). https://doi.org/10.1137/15M1054183; Shapeev, A. V. "Moment tensor potentials." *Multiscale Model. Simul.* **14**, 1153-1173 (2016). https://doi.org/10.1137/15M1054183
- Unified name-overlap candidate(s): MLIP ecosystem
- Unified closest-by-similarity (normalized): Molpro, i-PI, AMP

### MOLCAS
- IDs in entries tables: 91
- Context (nearest headings in source): Category: 1.1 PLANE-WAVE PSEUDOPOTENTIAL & PAW METHODS > Category: 1.2 ALL-ELECTRON & FULL-POTENTIAL METHODS > Category: 1.3 QUANTUM CHEMISTRY (GAUSSIAN BASIS)
- Source rows: `scientific_tools_consolidated/new_claude/entries_001_100_complete.md:120`
- License: Academic
- Official website: Predecessor to OpenMolcas
- Basis / method: Multireference
- Specialization: Legacy multireference package
- Primary publication (as in entries table): Historical predecessor to OpenMolcas, now superseded
- Unified name-overlap candidate(s): OpenMolcas
- Unified closest-by-similarity (normalized): OpenMolcas, LMTO-ASA, ORCA, JMol, TBmodels

### molecularGSM
- IDs in entries tables: 176, 234
- Context (nearest headings in source): Category: 5. PHONONS & LATTICE DYNAMICS > Category: 6. MOLECULAR DYNAMICS > Category: 7. TRANSITION STATES & RARE EVENTS; Sections 7.3 (Transition State Methods) & 8.1 (Structure Prediction) > Verified with 100% Accuracy Standard | January 2026 > Category: 7.3 Rare Events, Transition States & Enhanced Sampling
- Source rows: `scientific_tools_consolidated/new_claude/entries_101_200_complete.md:133`, `scientific_tools_consolidated/new_claude/entries_201_260_complete.md:416`, `scientific_tools_consolidated/new_claude/entries_234_237_transition_structure_prediction.md:9`
- License: MIT License
- Official website: https://github.com/ZimmermanGroup/molecularGSM; https://github.com/ZimmermanGroup/molecularGSM<br>https://zimmermangroup.github.io/molecularGSM/
- Basis / method: Growing String Method; Growing String Method for reaction paths
- Specialization: Double-ended (DE-GSM) and single-ended (SE-GSM) methods for transition state and reaction path finding, interfaces Gaussian/ORCA/TeraChem/Molpro/Q-Chem/MOPAC/ASE/xTB; Reaction pathway discovery, GSM
- Primary publication (as in entries table): Zimmerman, P. M. "Automated discovery of chemically reasonable elementary reaction steps." *J. Comput. Chem.* **34**, 1385-1392 (2013). https://doi.org/10.1002/jcc.23271; Zimmerman, P. M. "Reliable Transition State Searches Integrated with the Growing String Method." *J. Chem. Theory Comput.* **9** (7), 3043-3050 (2013). https://doi.org/10.1021/ct400319w<br><br>Zimmerman, P. M. "Single-Ended Transition State Finding with the Growing String Method." *J. Comput. Chem.* **36** (9), 601-611 (2015). https://doi.org/10.1002/jcc.23833

### NAMD
- IDs in entries tables: 170
- Context (nearest headings in source): Category: 4. TIGHT-BINDING & WANNIER FUNCTIONS > Category: 5. PHONONS & LATTICE DYNAMICS > Category: 6. MOLECULAR DYNAMICS
- Source rows: `scientific_tools_consolidated/new_claude/entries_101_200_complete.md:122`
- License: Free (academic)
- Official website: https://www.ks.uiuc.edu/Research/namd/
- Basis / method: Classical MD
- Specialization: NAnoscale Molecular Dynamics, biomolecular
- Primary publication (as in entries table): Phillips, J. C.; et al. "Scalable molecular dynamics on CPU and GPU architectures." *J. Chem. Phys.* **153**, 044130 (2020). https://doi.org/10.1063/5.0014475
- Unified closest-by-similarity (normalized): NOMAD, VMD, AMP, ALM, ADF

### NEGF-DFT
- IDs in entries tables: 117
- Context (nearest headings in source): Category: 1.5 ALL-ELECTRON NUMERIC BASIS > Category: 1.6 LINEAR-SCALING & ORDER-N METHODS > Category: 1.7 SPECIALIZED DFT IMPLEMENTATIONS
- Source rows: `scientific_tools_consolidated/new_claude/entries_101_200_complete.md:39`
- License: Research
- Official website: Various implementations
- Basis / method: Non-equilibrium Green's function
- Specialization: Transport, NEGF formalism
- Primary publication (as in entries table): Brandbyge, M.; et al. "Density-functional method for nonequilibrium electron transport." *Phys. Rev. B* **65**, 165401 (2002). https://doi.org/10.1103/PhysRevB.65.165401
- Unified closest-by-similarity (normalized): RMGDFT, EDMFTF, BigDFT

### NRG Ljubljana
- IDs in entries tables: 130
- Context (nearest headings in source): Category: 1.7 SPECIALIZED DFT IMPLEMENTATIONS > Category: 2. EXCITED-STATE METHODS (GW/BSE/TDDFT) > Category: 3. STRONGLY CORRELATED SYSTEMS
- Source rows: `scientific_tools_consolidated/new_claude/entries_101_200_complete.md:62`
- License: GNU GPL
- Official website: http://nrgljubljana.ijs.si/
- Basis / method: NRG
- Specialization: Numerical Renormalization Group, quantum impurity
- Primary publication (as in entries table): Žitko, R.; Pruschke, T. "Energy resolution and discretization artifacts in the numerical renormalization group." *Phys. Rev. B* **79**, 085106 (2009). https://doi.org/10.1103/PhysRevB.79.085106

### NumG
- IDs in entries tables: 107
- Context (nearest headings in source): Triple-Verified with 100% Accuracy Standard | January 2026 > Category: 1.4 GAUSSIAN BASIS - SPECIALIZED PACKAGES > Category: 1.5 ALL-ELECTRON NUMERIC BASIS
- Source rows: `scientific_tools_consolidated/new_claude/entries_101_200_complete.md:19`
- License: Research
- Official website: Limited distribution
- Basis / method: Numeric Gaussian
- Specialization: Numeric Gaussian approaches
- Primary publication (as in entries table): Research implementation of numeric Gaussian methods
- Unified closest-by-similarity (normalized): RMG

### OPTIMADE
- IDs in entries tables: 275
- Context (nearest headings in source): Triple-Verified with 100% Accuracy Standard | January 2026 > Category: 10.2 HIGH-THROUGHPUT & DATABASE INFRASTRUCTURE > Category: 10.2.2 SPECIALIZED HIGH-THROUGHPUT TOOLS
- Source rows: `scientific_tools_consolidated/new_claude/entries_261_290_complete.md:27`
- License: MIT License
- Official website: https://www.optimade.org/
- Basis / method: Materials database API
- Specialization: Open Databases Integration, unified API, 20+ databases
- Primary publication (as in entries table): Andersen, C. W.; Armiento, R.; Blokhin, E.; Conduit, G. J.; Dwaraknath, S.; Evans, M. L.; Fekete, Á.; Gopakumar, A.; et al. "OPTIMADE, an API for exchanging materials data." *Sci. Data* **8**, 217 (2021). https://doi.org/10.1038/s41597-021-00974-z
- Unified closest-by-similarity (normalized): NOMAD

### Pomerol
- IDs in entries tables: 214
- Context (nearest headings in source): Verified with 100% Accuracy Standard | January 2026 > Category: 4. WAVEFUNCTION-BASED QUANTUM CHEMISTRY > Category: 5. TIGHT-BINDING, MODEL HAMILTONIANS & DOWNFOLDING
- Source rows: `scientific_tools_consolidated/new_claude/entries_201_230_materials_codes.md:26`, `scientific_tools_consolidated/new_claude/entries_201_260_complete.md:26`
- License: GNU GPL v2
- Official website: https://github.com/aeantipov/pomerol
- Basis / method: Exact diagonalization solver
- Specialization: ED impurity solver, Green's functions, DMFT
- Primary publication (as in entries table): Antipov, A. E.; Krivenko, I.; Iskakov, S.; Kin-Lic Chan, G.; Rubtsov, A. N. "Pomerol: an exact diagonalization solver for the quantum impurity problem." (2014). https://github.com/aeantipov/pomerol (Code with documentation and examples)
- Unified name-overlap candidate(s): pomerol
- Unified closest-by-similarity (normalized): pomerol, PyMOL

### Priroda
- IDs in entries tables: 94
- Context (nearest headings in source): Category: 1.1 PLANE-WAVE PSEUDOPOTENTIAL & PAW METHODS > Category: 1.2 ALL-ELECTRON & FULL-POTENTIAL METHODS > Category: 1.3 QUANTUM CHEMISTRY (GAUSSIAN BASIS)
- Source rows: `scientific_tools_consolidated/new_claude/entries_001_100_complete.md:123`
- License: Academic
- Official website: Limited distribution
- Basis / method: QC
- Specialization: Russian QC development
- Primary publication (as in entries table): Laikov, D. N. "Fast evaluation of density functional theory." *Chem. Phys. Lett.* **281**, 151-156 (1997). https://doi.org/10.1016/S0009-2614(97)01206-2
- Unified name-overlap candidate(s): Priroda-06
- Unified closest-by-similarity (normalized): Priroda-06

### PySCF-DMRG
- IDs in entries tables: 140
- Context (nearest headings in source): Category: 2. EXCITED-STATE METHODS (GW/BSE/TDDFT) > Category: 3. STRONGLY CORRELATED SYSTEMS > Category: 3.2 MODEL HAMILTONIANS & PEDAGOGICAL
- Source rows: `scientific_tools_consolidated/new_claude/entries_101_200_complete.md:77`
- License: Part of PySCF
- Official website: DMRG-SCF
- Basis / method: Sharma, S.; Chan, G. K.-L. "Spin-adapted DMRG." *J. Chem. Phys.* **136**, 124121 (2012). https://doi.org/10.1063/1.3695642
- Primary publication (as in entries table): DMRG-SCF in PySCF
- Unified name-overlap candidate(s): PySCF
- Unified closest-by-similarity (normalized): PySCF

### PyXtal
- IDs in entries tables: 187
- Context (nearest headings in source): Category: 6. MOLECULAR DYNAMICS > Category: 7. TRANSITION STATES & RARE EVENTS > Category: 8. STRUCTURE PREDICTION
- Source rows: `scientific_tools_consolidated/new_claude/entries_101_200_complete.md:149`
- License: MIT License
- Official website: https://github.com/qzhu2017/PyXtal
- Basis / method: Crystal generator
- Specialization: Random crystal structure generation
- Primary publication (as in entries table): Fredericks, S.; et al. "PyXtal: A Python library for crystal structure generation." *Comput. Phys. Commun.* **261**, 107810 (2021). https://doi.org/10.1016/j.cpc.2020.107810
- Unified closest-by-similarity (normalized): XtalOpt, CRYSTAL

### QuSpin
- IDs in entries tables: 138
- Context (nearest headings in source): Category: 2. EXCITED-STATE METHODS (GW/BSE/TDDFT) > Category: 3. STRONGLY CORRELATED SYSTEMS > Category: 3.2 MODEL HAMILTONIANS & PEDAGOGICAL
- Source rows: `scientific_tools_consolidated/new_claude/entries_101_200_complete.md:75`
- License: BSD 3-Clause
- Official website: https://weinbe58.github.io/QuSpin/
- Basis / method: Exact diagonalization
- Specialization: Python ED, spin systems, dynamics
- Primary publication (as in entries table): Weinberg, P.; Bukov, M. "QuSpin: a Python package for dynamics and exact diagonalization." *SciPost Phys.* **2**, 003 (2017). https://doi.org/10.21468/SciPostPhys.2.1.003
- Unified closest-by-similarity (normalized): S/PHI/nX, Gaussian

### Random Search
- IDs in entries tables: 190
- Context (nearest headings in source): Category: 6. MOLECULAR DYNAMICS > Category: 7. TRANSITION STATES & RARE EVENTS > Category: 8. STRUCTURE PREDICTION
- Source rows: `scientific_tools_consolidated/new_claude/entries_101_200_complete.md:152`
- License: Various
- Official website: Various implementations
- Basis / method: Random sampling
- Specialization: Random structure sampling
- Primary publication (as in entries table): Various random structure generation implementations

### SlateKoster
- IDs in entries tables: 149
- Context (nearest headings in source): Category: 3. STRONGLY CORRELATED SYSTEMS > Category: 3.2 MODEL HAMILTONIANS & PEDAGOGICAL > Category: 4. TIGHT-BINDING & WANNIER FUNCTIONS
- Source rows: `scientific_tools_consolidated/new_claude/entries_101_200_complete.md:91`
- License: Research
- Official website: Various implementations
- Basis / method: Slater-Koster
- Specialization: Slater-Koster parameterization
- Primary publication (as in entries table): Slater, J. C.; Koster, G. F. "Simplified LCAO method." *Phys. Rev.* **94**, 1498-1524 (1954). https://doi.org/10.1103/PhysRev.94.1498
- Unified closest-by-similarity (normalized): Lobster

### SMEAGOL
- IDs in entries tables: 116
- Context (nearest headings in source): Category: 1.5 ALL-ELECTRON NUMERIC BASIS > Category: 1.6 LINEAR-SCALING & ORDER-N METHODS > Category: 1.7 SPECIALIZED DFT IMPLEMENTATIONS
- Source rows: `scientific_tools_consolidated/new_claude/entries_101_200_complete.md:38`
- License: Academic
- Official website: https://www.smeagol.tcd.ie/
- Basis / method: DFT + transport
- Specialization: Spin and Molecular Electronics in Atomically Generated Orbital Landscapes
- Primary publication (as in entries table): Rocha, A. R.; et al. "Towards molecular spintronics." *Nat. Mater.* **4**, 335-339 (2005). https://doi.org/10.1038/nmat1349
- Unified closest-by-similarity (normalized): pomerol

### SNAP
- IDs in entries tables: 319
- Context (nearest headings in source): Category: 11.6 SPECIALIZED SOLVERS & METHODS > Category: 11.7 ADDITIONAL SPECIALIZED CODES > Category: 11.7.2 ADDITIONAL FRAMEWORKS & UTILITIES
- Source rows: `scientific_tools_consolidated/new_claude/entries_291_320_complete.md:46`
- License: BSD-like
- Official website: https://github.com/FitSNAP/FitSNAP
- Basis / method: ML potential fitting
- Specialization: Spectral neighbor analysis, bispectrum descriptors, linear ML potentials, Sandia
- Primary publication (as in entries table): Thompson, A. P.; Swiler, L. P.; Trott, C. R.; Foiles, S. M.; Tucker, G. J. "Spectral neighbor analysis method for automated generation of quantum-accurate interatomic potentials." *J. Comput. Phys.* **285**, 316-330 (2015). https://doi.org/10.1016/j.jcp.2014.12.018; Wood, M. A.; Thompson, A. P. "Extending the accuracy of the SNAP interatomic potential form." *J. Chem. Phys.* **148**, 241721 (2018). https://doi.org/10.1063/1.5017641
- Unified closest-by-similarity (normalized): AMP

### SPHInX
- IDs in entries tables: 98
- Context (nearest headings in source): Category: 1.1 PLANE-WAVE PSEUDOPOTENTIAL & PAW METHODS > Category: 1.2 ALL-ELECTRON & FULL-POTENTIAL METHODS > Category: 1.3 QUANTUM CHEMISTRY (GAUSSIAN BASIS)
- Source rows: `scientific_tools_consolidated/new_claude/entries_001_100_complete.md:127`
- License: Academic
- Official website: https://sxrepo.mpie.de/
- Basis / method: Numeric basis DFT
- Specialization: S/PHI/nX numeric AO, Max Planck
- Primary publication (as in entries table): Boeck, S.; et al. "The object-oriented DFT program library S/PHI/nX." *Comput. Phys. Commun.* **182**, 543-554 (2011). https://doi.org/10.1016/j.cpc.2010.09.016
- Unified closest-by-similarity (normalized): S/PHI/nX, Spex, PHON

### SPINPACK
- IDs in entries tables: 137
- Context (nearest headings in source): Category: 2. EXCITED-STATE METHODS (GW/BSE/TDDFT) > Category: 3. STRONGLY CORRELATED SYSTEMS > Category: 3.2 MODEL HAMILTONIANS & PEDAGOGICAL
- Source rows: `scientific_tools_consolidated/new_claude/entries_101_200_complete.md:74`
- License: Research
- Official website: Various implementations
- Basis / method: Spin models
- Specialization: Heisenberg, Ising, XY models
- Primary publication (as in entries table): Various implementations of spin model solvers
- Unified closest-by-similarity (normalized): SchNetPack, EDIpack, Z2Pack, S/PHI/nX

### StoBe
- IDs in entries tables: 115
- Context (nearest headings in source): Category: 1.5 ALL-ELECTRON NUMERIC BASIS > Category: 1.6 LINEAR-SCALING & ORDER-N METHODS > Category: 1.7 SPECIALIZED DFT IMPLEMENTATIONS
- Source rows: `scientific_tools_consolidated/new_claude/entries_101_200_complete.md:37`
- License: Academic
- Official website: http://www.fhi-berlin.mpg.de/KHsoftware/StoBe/
- Basis / method: STOs
- Specialization: Stockholm-Berlin DFT, STOs
- Primary publication (as in entries table): Herrmann, C.; et al. "StoBe-deMon." (Code documentation, FHI Berlin). http://www.fhi-berlin.mpg.de/KHsoftware/StoBe/
- Unified closest-by-similarity (normalized): Stoner

### String Method
- IDs in entries tables: 179
- Context (nearest headings in source): Category: 5. PHONONS & LATTICE DYNAMICS > Category: 6. MOLECULAR DYNAMICS > Category: 7. TRANSITION STATES & RARE EVENTS
- Source rows: `scientific_tools_consolidated/new_claude/entries_101_200_complete.md:136`
- License: Research
- Official website: Various implementations
- Basis / method: String methods
- Specialization: Minimum energy pathways
- Primary publication (as in entries table): E, W.; et al. "String method for the study of rare events." *Phys. Rev. B* **66**, 052301 (2002). https://doi.org/10.1103/PhysRevB.66.052301
- Unified name-overlap candidate(s): String methods
- Unified closest-by-similarity (normalized): String methods

### TBE
- IDs in entries tables: 148
- Context (nearest headings in source): Category: 3. STRONGLY CORRELATED SYSTEMS > Category: 3.2 MODEL HAMILTONIANS & PEDAGOGICAL > Category: 4. TIGHT-BINDING & WANNIER FUNCTIONS
- Source rows: `scientific_tools_consolidated/new_claude/entries_101_200_complete.md:90`
- License: Research
- Official website: Various implementations
- Basis / method: Tight-binding
- Specialization: General TB codes
- Primary publication (as in entries table): Various tight-binding eigenvalue solvers
- Unified name-overlap candidate(s): Matbench
- Unified closest-by-similarity (normalized): xTB, TDEP, TB2J, NBSE, KITE

### TenPy
- IDs in entries tables: 132
- Context (nearest headings in source): Category: 2. EXCITED-STATE METHODS (GW/BSE/TDDFT) > Category: 3. STRONGLY CORRELATED SYSTEMS > Category: 3.2 MODEL HAMILTONIANS & PEDAGOGICAL
- Source rows: `scientific_tools_consolidated/new_claude/entries_101_200_complete.md:69`
- License: GNU GPL v3
- Official website: https://tenpy.readthedocs.io/
- Basis / method: Tensor networks
- Specialization: Python tensor networks, DMRG, MPS
- Primary publication (as in entries table): Hauschild, J.; Pollmann, F. "Efficient numerical simulations with Tensor Networks." *SciPost Phys. Lect. Notes* **5** (2018). https://doi.org/10.21468/SciPostPhysLectNotes.5
- Unified closest-by-similarity (normalized): TDEP, MatPy

### TINKER
- IDs in entries tables: 175
- Context (nearest headings in source): Category: 4. TIGHT-BINDING & WANNIER FUNCTIONS > Category: 5. PHONONS & LATTICE DYNAMICS > Category: 6. MOLECULAR DYNAMICS
- Source rows: `scientific_tools_consolidated/new_claude/entries_101_200_complete.md:127`
- License: Custom (academic)
- Official website: https://dasher.wustl.edu/tinker/
- Basis / method: Classical MD
- Specialization: Molecular design, force field development
- Primary publication (as in entries table): Ponder, J. W. "TINKER: Software Tools for Molecular Design." (Washington University 2004+). https://dasher.wustl.edu/tinker/
- Unified closest-by-similarity (normalized): matminer, Stoner

### TranSIESTA
- IDs in entries tables: 118
- Context (nearest headings in source): Category: 1.5 ALL-ELECTRON NUMERIC BASIS > Category: 1.6 LINEAR-SCALING & ORDER-N METHODS > Category: 1.7 SPECIALIZED DFT IMPLEMENTATIONS
- Source rows: `scientific_tools_consolidated/new_claude/entries_101_200_complete.md:40`
- License: GNU GPL
- Official website: Part of SIESTA
- Basis / method: DFT + transport
- Specialization: NEGF transport with SIESTA
- Primary publication (as in entries table): Brandbyge, M.; et al. "Density-functional method for nonequilibrium electron transport." *Phys. Rev. B* **65**, 165401 (2002). https://doi.org/10.1103/PhysRevB.65.165401
- Unified name-overlap candidate(s): SIESTA
- Unified closest-by-similarity (normalized): SIESTA, fiesta

### USPEX-ML
- IDs in entries tables: 188
- Context (nearest headings in source): Category: 6. MOLECULAR DYNAMICS > Category: 7. TRANSITION STATES & RARE EVENTS > Category: 8. STRUCTURE PREDICTION
- Source rows: `scientific_tools_consolidated/new_claude/entries_101_200_complete.md:150`
- License: Academic
- Official website: Part of USPEX
- Basis / method: ML-enhanced CSP
- Specialization: ML-accelerated structure prediction
- Primary publication (as in entries table): Podryabinkin, E. V.; et al. "Accelerating crystal structure prediction by machine-learning interatomic potentials." *Phys. Rev. B* **99**, 064114 (2019). https://doi.org/10.1103/PhysRevB.99.064114
- Unified name-overlap candidate(s): Spex, USPEX
- Unified closest-by-similarity (normalized): USPEX, Spex

### VASP-BSE
- IDs in entries tables: 124
- Context (nearest headings in source): Category: 1.6 LINEAR-SCALING & ORDER-N METHODS > Category: 1.7 SPECIALIZED DFT IMPLEMENTATIONS > Category: 2. EXCITED-STATE METHODS (GW/BSE/TDDFT)
- Source rows: `scientific_tools_consolidated/new_claude/entries_101_200_complete.md:51`
- License: Commercial
- Official website: Part of VASP
- Basis / method: BSE
- Specialization: BSE in VASP
- Primary publication (as in entries table): Shishkin, M.; Kresse, G. "Implementation of the GW method for semiconductors." *Phys. Rev. B* **74**, 035101 (2006). https://doi.org/10.1103/PhysRevB.74.035101
- Unified name-overlap candidate(s): VASP
- Unified closest-by-similarity (normalized): VASP, ASE, vaspkit, VAMPIRE

### VASP-GW
- IDs in entries tables: 49
- Context (nearest headings in source): Category: 1.1 PLANE-WAVE PSEUDOPOTENTIAL & PAW METHODS > Category: 1.2 ALL-ELECTRON & FULL-POTENTIAL METHODS > Category: 1.3 QUANTUM CHEMISTRY (GAUSSIAN BASIS)
- Source rows: `scientific_tools_consolidated/new_claude/entries_001_100_complete.md:69`
- License: Commercial
- Official website: Part of VASP
- Basis / method: GW
- Specialization: GW implementation in VASP
- Primary publication (as in entries table): Shishkin, M.; Kresse, G. "Implementation and performance of the frequency-dependent GW method." *Phys. Rev. B* **74**, 035101 (2006). https://doi.org/10.1103/PhysRevB.74.035101
- Unified name-overlap candidate(s): VASP
- Unified closest-by-similarity (normalized): VASP, vaspkit, GASP, VASP+DMFT

## Extra In Unified (not in entries*.md)

- ACES II (overlap: ACE, ACES; similar: ACES, ACE, MACE)
- aiida-fleur (overlap: AiiDA, Fleur; similar: Fleur, AiiDA, AiiDA-wannier90)
- ALF (similar: ALM, ADF, ALPS)
- ALPS/CT-HYB (overlap: ALPS; similar: ALPS)
- ALPSCore (overlap: ALPS; similar: CALYPSO, ALPS, DCore)
- AMULET (similar: AMSET)
- Atomistix ToolKit
- Atompaw/PWPAW (overlap: Atompaw, AtomPAW; similar: Atompaw, AtomPAW)
- AtomViz (similar: FSatom, atomate, Atompaw, AtomPAW)
- Basin hopping (similar: Pybinding, CASINO)
- CHAMP (overlap: AMP; similar: AMP, CatMAP, CHARMM, cmpy, COHP)
- CT-HYB
- CT-INT (similar: exciting, custodian)
- CT-QMC
- CT-SEG
- DCA++ (similar: ADC, ORCA, DDEC, CCDC)
- DMDW/RTDW (overlap: DMDW; similar: DMDW)
- DMFTwDFT (similar: EDMFTF, RMGDFT)
- EPW (similar: TDEP, GPAW)
- fiesta (similar: SIESTA, VESTA, TranSIESTA, WEST, Questaal)
- FTPS (similar: PQS)
- GAMESS-UK (similar: GAMESS-US)
- Gaussian 16 (overlap: Gaussian; similar: Gaussian)
- Gaussian Basis Set Libraries (overlap: Gaussian)
- GreenX (similar: FreeON)
- HOTBIT (similar: HTST)
- HÏ†
- jobflow-remote (overlap: jobflow; similar: jobflow)
- JuKKR
- Julia materials packages (similar: Julia packages, Materials Project, Materials Cloud)
- KKRnano
- LMTO-ASA
- Materials Simulation Suites (similar: Materials and Processes Simulations, Materials Studio)
- Metadynamics (similar: w2dynamics)
- MLIP ecosystem (overlap: MLIP)
- molgw (similar: JMol, PyMOL, DMOL3)
- n2p2
- NBSE (overlap: BSE; similar: BSE, TBE, NEB, ASE)
- Neural network potentials
- PLATO (similar: Petot, gpaw-tools)
- Priroda-06 (overlap: Priroda; similar: Priroda)
- pymatgen-diffusion (overlap: pymatgen; similar: pymatgen-db, pymatgen)
- PyQMC (similar: MPQC, PySCF, PyMOL)
- PyQuante
- Qbox
- QMcBeaver (similar: Bader, Amber)
- QUEST (overlap: CONQUEST, Questaal; similar: Questaal, CONQUEST, WEST, iQIST, VESTA)
- QWalk
- RMG (overlap: RMGDFT; similar: RMGDFT, NumG, NORG, DMRG++)
- SALMON (overlap: ALM; similar: Dalton, ALM, ALAMODE, sumo, dual fermions)
- SchNetPack (similar: SPINPACK, CHGNet, QMCPACK, EDIpack)
- solid_dmft (similar: ComDMFT)
- Spex (overlap: USPEX, USPEX-ML; similar: USPEX, USPEX-ML, COPEX, SPHInX, S/PHI/nX)
- Stoner (similar: StoBe, TINKER, Lobster, matminer)
- String methods (overlap: String Method; similar: String Method)
- TOTAL (similar: Petot)
- TRIQS/CT-QMC solvers (overlap: TRIQS)
- TRIQS/DFTTools (overlap: TRIQS; similar: TRIQS)
- VASP+DMFT (overlap: VASP; similar: vaspkit, VASP, VASP-GW, EDMFTF)
- X-ray Crystallography Integration (overlap: CRYSTAL)
- Z2Pack (similar: QMCPACK, EDIpack, SPINPACK)
- ~141
- ~328
- ~469
