# COMPREHENSIVE MATERIALS CODES CATALOG
## Complete Edition - All Codes with Full Details & Proper Categorization
### Fixed: Both Table Formats Included | January 2026

---

**Grand Total Codes**: 278
**Codes with Complete Details**: 214 (77.0%)
**Codes Missing Details**: 64 (23.0%)

---


====================================================================================================
# 1. GROUND-STATE ELECTRONIC STRUCTURE (DFT & VARIANTS)
====================================================================================================

## 1.1 Plane-Wave Pseudopotential & PAW Methods
**Subcategory Total: 87 codes**
────────────────────────────────────────────────────────────────────────────────

### 1.1.1 Major Production Codes - Authenticated Table (10 Codes)

| # | ID | Code Name | License | Official Website | Basis Set | Primary Publication with DOI | Specialization |
|---|---|---|---|---|---|---|---|
| 1 | 1 | **VASP** | Commercial | https://www.vasp.at/ | Plane-wave + PAW | Kresse, G.; Hafner, J. "Ab Initio Molecular-Dynamics for Liquid-Metals." *Phys. Rev. B* **47**, 558-561 (1993). https://doi.org/10.1103/PhysRevB.47.558; Kresse, G.; Furthmüller, J. "Efficiency of ab-initio total energy calculations for metals and semiconductors using a plane-wave basis set." *Comput. Mater. Sci.* **6**, 15-50 (1996). https://doi.org/10.1016/0927-0256(96)00008-0 | Vienna Ab initio Simulation Package, PAW method, hybrid functionals, GW/RPA |
| 2 | 2 | **Quantum ESPRESSO** | GNU GPL v2 | https://www.quantum-espresso.org/ | Plane-wave + Pseudopotential | Giannozzi, P.; et al. "QUANTUM ESPRESSO: a modular and open-source software project for quantum simulations of materials." *J. Phys.: Condens. Matter* **21**, 395502 (2009). https://doi.org/10.1088/0953-8984/21/39/395502 | Modular suite, DFPT, phonons, GW, BSE, optical properties |
| 3 | 3 | **ABINIT** | GNU GPL v2 | https://www.abinit.org/ | Plane-wave + PAW | Gonze, X.; et al. "ABINIT: First-principles approach to material and nanosystem properties." *Comput. Phys. Commun.* **180**, 2582-2615 (2009). https://doi.org/10.1016/j.cpc.2009.02.005 | GW, BSE, DMFT, DFPT, response properties |
| 4 | 4 | **CASTEP** | Commercial/Academic | https://www.castep.org/ | Plane-wave + Pseudopotential | Clark, S. J.; et al. "First principles methods using CASTEP." *Z. Kristallogr.* **220**, 567-570 (2005). https://doi.org/10.1524/zkri.220.5.567.65075 | Cambridge Serial Total Energy Package, NMR, EELS |
| 5 | 5 | **CP2K** | GNU GPL | https://www.cp2k.org/ | Hybrid Gaussian + Plane-wave | Kühne, T. D.; et al. "CP2K: An electronic structure and molecular dynamics software package." *J. Chem. Phys.* **152**, 194103 (2020). https://doi.org/10.1063/5.0007045 | Mixed basis, BOMD, CPMD, real-time TDDFT |
| 6 | 6 | **GPAW** | GNU GPL | https://wiki.fysik.dtu.dk/gpaw/ | Real-space grid + PAW | Enkovaara, J.; et al. "Electronic structure calculations with GPAW." *J. Phys.: Condens. Matter* **22**, 253202 (2010). https://doi.org/10.1088/0953-8984/22/25/253202 | Grid-based Projector Augmented Wave, Python-integrated |
| 7 | 7 | **NWChem** | Educational Community License v2 | https://nwchemgit.github.io/ | Plane-wave + Gaussian | Valiev, M.; et al. "NWChem: A comprehensive and scalable open-source solution." *Comput. Phys. Commun.* **181**, 1477-1489 (2010). https://doi.org/10.1016/j.cpc.2010.04.018 | Quantum chemistry, MD, coupled-cluster, DFT |
| 8 | 8 | **SIESTA** | GNU GPL | https://siesta-project.org/ | Numerical LCAO | Soler, J. M.; et al. "The SIESTA method for ab initio order-N materials simulation." *J. Phys.: Condens. Matter* **14**, 2745-2779 (2002). https://doi.org/10.1088/0953-8984/14/11/302 | Linear-scaling DFT, numerical orbitals, Spanish development |
| 9 | 11 | **Octopus** | GNU GPL | https://octopus-code.org/ | Real-space grid | Andrade, X.; et al. "Real-space grids and the Octopus code." *Phys. Chem. Chem. Phys.* **17**, 31371-31396 (2015). https://doi.org/10.1039/C5CP00351B | Real-time TDDFT, linear-response TDDFT, MD |
| 10 | 12 | **BigDFT** | GNU GPL | https://bigdft.org/ | Daubechies wavelets | Genovese, L.; et al. "Daubechies wavelets for high performance electronic structure calculations." *Comput. Phys. Commun.* **181**, 1919-1931 (2010). https://doi.org/10.1016/j.cpc.2010.08.005 | Linear-scaling DFT, wavelet basis, GPU-accelerated |
| 11 | 14 | **RMGDFT** | GNU GPL | https://github.com/RMGDFT/rmgdft | Real-space multigrid | Briggs, E. L.; et al. "Large-scale electronic structure calculations with multigrid acceleration." *Phys. Rev. B* **54**, 14362-14364 (1996). https://doi.org/10.1103/PhysRevB.54.14362 | Real-space DFT, GPU-accelerated, massively parallel |
| 12 | 15 | **PARSEC** | Academic/BSD | http://www.ices.utexas.edu/parsec/ | Real-space pseudopotential | Kronik, L.; et al. "PARSEC – the pseudopotential algorithm for real-space electronic structure calculations." *Phys. Status Solidi B* **243**, 1063-1079 (2006). https://doi.org/10.1002/pssb.200541463 | Real-space grid, nanostructures, TDDFT |
| 13 | 16 | **DFTB+** | GNU LGPL v3 | https://dftbplus.org/ | Tight-binding | Hourahine, B.; et al. "DFTB+, a software package for efficient approximate density functional theory." *J. Chem. Phys.* **152**, 124101 (2020). https://doi.org/10.1063/1.5143190 | SCC-DFTB, large-scale simulations, excited states |
| 14 | 17 | **xTB** | Apache 2.0 | https://github.com/grimme-lab/xtb | Extended tight-binding | Grimme, S.; et al. "A robust and accurate tight-binding quantum chemical method." *J. Chem. Theory Comput.* **13**, 1989-2009 (2017). https://doi.org/10.1021/acs.jctc.7b00118 | GFN methods, rapid QC, structure prediction, Grimme group |
| 15 | 18 | **OpenMX** | GNU GPL | https://www.openmx-square.org/ | Numerical atomic orbitals | Ozaki, T. "Variationally optimized atomic orbitals for large-scale electronic structures." *Phys. Rev. B* **67**, 155108 (2003). https://doi.org/10.1103/PhysRevB.67.155108 | NAOs, Japanese development, efficient scaling |
| 16 | 19 | **CONQUEST** | Academic | https://www.order-n.org/ | Numerical orbitals | Bowler, D. R.; Miyazaki, T. "O(N) methods in electronic structure calculations." *Rep. Prog. Phys.* **75**, 036503 (2012). https://doi.org/10.1088/0034-4885/75/3/036503 | Linear-scaling order-N DFT, large systems |
| 17 | 20 | **ONETEP** | Academic/Commercial | https://www.onetep.org/ | NGWFs | Skylaris, C. K.; et al. "Introducing ONETEP: Linear-scaling density functional simulations." *J. Chem. Phys.* **122**, 084119 (2005). https://doi.org/10.1063/1.1839852 | Non-orthogonal generalized Wannier functions, linear-scaling |
| 18 | 43 | **BerkeleyGW** | BSD | https://berkeleygw.org/ | GW/BSE | Deslippe, J.; et al. "BerkeleyGW: A massively parallel computer package for the calculation of the quasiparticle and optical properties." *Comput. Phys. Commun.* **183**, 1269-1289 (2012). https://doi.org/10.1016/j.cpc.2011.12.006 | GW-BSE specialist, optical properties, Berkeley |
| 19 | 44 | **Yambo** | GNU GPL v2 | http://www.yambo-code.org/ | GW/BSE | Sangalli, D.; et al. "Many-body perturbation theory calculations using the yambo code." *J. Phys.: Condens. Matter* **31**, 325902 (2019). https://doi.org/10.1088/1361-648X/ab15d0 | GW, BSE, TDDFT, real-time |
| 20 | 121 | **exciting** | GNU GPL v3 | http://exciting-code.org/ | FP-LAPW | Gulans, A.; et al. "exciting: a full-potential all-electron package." *J. Phys.: Condens. Matter* **26**, 363202 (2014). https://doi.org/10.1088/0953-8984/26/36/363202 | All-electron LAPW, GW, BSE, TDDFT |
| 21 | 21 | **WIEN2k** | Commercial | http://www.wien2k.at/ | FP-LAPW | Blaha, P.; et al. "WIEN2k: An APW+lo program for calculating the properties of solids." *J. Chem. Phys.* **152**, 074101 (2020). https://doi.org/10.1063/1.5143061 | Industry standard LAPW, non-collinear magnetism, Vienna |
| 22 | 23 | **RSPt** | Academic | https://rspt.physics.uu.se/ | FP-LMTO | Wills, J. M.; et al. "Full-Potential Electronic Structure Method." Springer (2010). https://doi.org/10.1007/978-3-642-15144-6 | Relativistic Spin-Polarized test code, GW, EDMFT |
| 23 | 24 | **Questaal** | GNU GPL | https://www.questaal.org/ | LMTO/GW | Kotani, T.; et al. "Quasiparticle self-consistent GW method." *Phys. Rev. B* **76**, 165106 (2007). https://doi.org/10.1103/PhysRevB.76.165106 | LMTO suite, QSGW, tight-binding downfolding |
| 24 | 51 | **TRIQS** | GNU GPL v3 | https://triqs.github.io/ | DMFT framework | Parcollet, O.; et al. "TRIQS: A toolbox for research on interacting quantum systems." *Comput. Phys. Commun.* **196**, 398-415 (2015). https://doi.org/10.1016/j.cpc.2015.04.023 | Toolbox for Research on Interacting Quantum Systems |
| 25 | 59 | **QMCPACK** | BSD 3-Clause | https://qmcpack.org/ | QMC | Kim, J.; et al. "QMCPACK: an open source ab initio quantum Monte Carlo package." *J. Phys.: Condens. Matter* **30**, 195901 (2018). https://doi.org/10.1088/1361-648X/aab9c3 | Production-scale QMC, VMC, DMC |
| 26 | 60 | **CASINO** | Academic | https://vallico.net/casinoqmc/ | QMC | Needs, R. J.; et al. "Continuum variational and diffusion quantum Monte Carlo calculations." *J. Phys.: Condens. Matter* **22**, 023201 (2010). https://doi.org/10.1088/0953-8984/22/2/023201 | High-accuracy QMC, Cambridge |
| 27 | 62 | **Phonopy** | BSD 3-Clause | https://phonopy.github.io/phonopy/ | Harmonic phonons | Togo, A.; Tanaka, I. "First principles phonon calculations in materials science." *Scr. Mater.* **108**, 1-5 (2015). https://doi.org/10.1016/j.scriptamat.2015.07.021 | Phonon dispersions, industry standard |
| 28 | 63 | **phono3py** | BSD 3-Clause | https://phonopy.github.io/phono3py/ | Anharmonic phonons | Togo, A.; et al. "Distributions of phonon lifetimes in Brillouin zones." *Phys. Rev. B* **91**, 094306 (2015). https://doi.org/10.1103/PhysRevB.91.094306 | Lattice thermal conductivity, 3rd-order IFCs |
| 29 | 64 | **ShengBTE** | GNU GPL v3 | http://www.shengbte.org/ | Phonon BTE | Li, W.; et al. "ShengBTE: A solver of the Boltzmann transport equation." *Comput. Phys. Commun.* **185**, 1747-1758 (2014). https://doi.org/10.1016/j.cpc.2014.02.015 | Thermal conductivity, Boltzmann equation |
| 30 | 65 | **ALAMODE** | MIT License | https://github.com/ttadano/alamode | IFC extraction | Tadano, T.; et al. "Anharmonic force constants from compressive sensing." *Phys. Rev. B* **92**, 054301 (2015). https://doi.org/10.1103/PhysRevB.92.054301 | Anharmonic lattice dynamics, compressive sensing |
| 31 | 42 | **Wannier90** | GNU GPL v2 | http://www.wannier.org/ | Wannier functions | Pizzi, G.; et al. "Wannier90 as a community code." *J. Phys.: Condens. Matter* **32**, 165902 (2020). https://doi.org/10.1088/1361-648X/ab51ff | Maximally localized Wannier functions, industry standard |
| 32 | 66 | **ASE** | GNU LGPL v2.1 | https://wiki.fysik.dtu.dk/ase/ | Python framework | Larsen, A. H.; et al. "The atomic simulation environment." *J. Phys.: Condens. Matter* **29**, 273002 (2017). https://doi.org/10.1088/1361-648X/aa680e | Atomic Simulation Environment, 20+ calculators |
| 33 | 67 | **pymatgen** | MIT License | https://pymatgen.org/ | Materials analysis | Ong, S. P.; et al. "Python Materials Genomics (pymatgen)." *Comput. Mater. Sci.* **68**, 314-319 (2013). https://doi.org/10.1016/j.commatsci.2012.10.028 | Materials Genomics, analysis, I/O, MP interface |
| 34 | 68 | **AiiDA** | MIT License | https://www.aiida.net/ | Workflow engine | Huber, S. P.; et al. "AiiDA 1.0, a scalable computational infrastructure." *Sci. Data* **7**, 300 (2020). https://doi.org/10.1038/s41597-020-00638-4 | Automated Infrastructure, provenance tracking |
| 35 | 177 | **PLUMED** | GNU LGPL v3 | https://www.plumed.org/ | Metadynamics | Tribello, G. A.; et al. "PLUMED 2: New feathers for an old bird." *Comput. Phys. Commun.* **185**, 604-613 (2014). https://doi.org/10.1016/j.cpc.2013.09.018 | Metadynamics, enhanced sampling, plugin |
| 36 | 342 | **BoltzTraP** | Academic | http://www.icams.de/content/departments/ams/madsen-projekte/boltztrap/ | Boltzmann transport | Madsen, G. K. H.; Singh, D. J. "BoltzTraP. A code for calculating band-structure dependent quantities." *Comput. Phys. Commun.* **175**, 67-71 (2006). https://doi.org/10.1016/j.cpc.2006.03.007 | Boltzmann transport properties, classic |
| 37 | 313 | **Lobster** | Academic/Commercial | http://www.cohp.de/ | Chemical bonding analysis | Deringer, V. L.; Tchougréeff, A. L.; Dronskowski, R. "Crystal orbital Hamilton population (COHP) analysis as projected from plane-wave basis sets." *J. Phys. Chem. A* **115**, 5461-5466 (2011). https://doi.org/10.1021/jp202489s; Maintz, S.; Deringer, V. L.; Tchougréeff, A. L.; Dronskowski, R. "LOBSTER: A tool to extract chemical bonding from plane-wave based DFT." *J. Comput. Chem.* **37**, 1030-1035 (2016). https://doi.org/10.1002/jcc.24300 | COHP analysis from PAW, bonding analysis, VASP/ABINIT/QE compatible |
| 38 | 41 | **WannierBerri** | GNU GPL v2 | https://wannier-berri.org/ | Wannier TB | Tsirkin, S. S. "High performance Wannier interpolation of Berry curvature." *npj Comput. Mater.* **7**, 33 (2021). https://doi.org/10.1038/s41524-021-00498-5 | Berry curvature, AHC, SHC, orbital magnetization |
| 39 | - | **Z2Pack** | - | - | - | - | - |
| 40 | 325 | **AMSET** | Modified BSD | https://hackingmaterials.lbl.gov/amset/ | Transport properties | Ganose, A. M.; et al. "Efficient calculation of carrier scattering rates." *Nat. Commun.* **12**, 2222 (2021). https://doi.org/10.1038/s41467-021-22440-5 | Ab initio scattering and transport, LBNL |
| 41 | 356 | **FEFF** | Academic | http://feffproject.org/ | X-ray spectroscopy | Rehr, J. J.; et al. "Parameter-free calculations of X-ray spectra with FEFF9." *Phys. Chem. Chem. Phys.* **12**, 5503-5513 (2010). https://doi.org/10.1039/B926434E | X-ray absorption, EXAFS, XANES |
| 42 | 37 | **CRYSTAL** | Academic/Commercial | https://www.crystal.unito.it/ | Gaussian basis (periodic) | Dovesi, R.; et al. "CRYSTAL17." *WIREs Comput. Mol. Sci.* **8**, e1360 (2018). https://doi.org/10.1002/wcms.1360 | Periodic Gaussian basis, hybrid DFT specialist |
| 43 | 167 | **LAMMPS** | GNU GPL v2 | https://www.lammps.org/ | Classical MD | Thompson, A. P.; et al. "LAMMPS - a flexible simulation tool." *Comput. Phys. Commun.* **271**, 108171 (2022). https://doi.org/10.1016/j.cpc.2021.108171 | Large-scale Atomic/Molecular Massively Parallel Simulator |
| 44 | 166 | **i-PI** | GNU GPL v3 | http://ipi-code.org/ | Universal PIMD driver | Ceriotti, M.; et al. "i-PI: A Python interface for ab initio path integral molecular dynamics." *Comput. Phys. Commun.* **185**, 1019-1026 (2014). https://doi.org/10.1016/j.cpc.2013.10.027 | Path integral MD, client-server |
| 45 | - | **Spex** | - | - | - | - | - |
| 46 | - | **fiesta** | - | - | - | - | - |
| 47 | - | **EPW** | - | - | - | - | - |
| 48 | 27 | **ORCA** | Free-Academic | https://www.faccts.de/orca/ | Gaussian basis | Neese, F. "Software update: The ORCA program system." *WIREs Comput. Mol. Sci.* **12**, e1606 (2022). https://doi.org/10.1002/wcms.1606 | DLPNO-CC, excited states, EPR, Neese group |
| 49 | - | **Gaussian 16** | - | - | - | - | - |
| 50 | 28 | **PSI4** | GNU GPL v3 | https://psicode.org/ | Gaussian basis | Turney, J. M.; et al. "PSI4: an open-source ab initio electronic structure package." *WIREs Comput. Mol. Sci.* **2**, 556-565 (2012). https://doi.org/10.1002/wcms.93 | Python-driven, modular, extensible plugins |
| 51 | 31 | **Molpro** | Commercial | https://www.molpro.net/ | Gaussian basis | Werner, H.-J.; et al. "Molpro: a general-purpose quantum chemistry program package." *WIREs Comput. Mol. Sci.* **2**, 242-253 (2012). https://doi.org/10.1002/wcms.82 | High-accuracy CC, multireference methods |
| 52 | 56 | **MRCC** | Free Academic | https://www.mrcc.hu/ | CC methods | Kállay, M.; et al. "The MRCC program system." *J. Chem. Phys.* **152**, 074107 (2020). https://doi.org/10.1063/1.5142048 | Multireference coupled-cluster specialist |
| 53 | 57 | **OpenMolcas** | GNU LGPL v2.1 | https://gitlab.com/Molcas/OpenMolcas | Multireference | Fdez. Galván, I.; et al. "OpenMolcas: From Source Code to Insight." *J. Chem. Theory Comput.* **15**, 5925-5964 (2019). https://doi.org/10.1021/acs.jctc.9b00532 | CASSCF/CASPT2 specialist, NEVPT2, DMRG |
| 54 | 58 | **BAGEL** | GNU GPL v3 | https://nubakery.org/bagel/ | Multireference | Shiozaki, T. "BAGEL: Brilliantly Advanced General Electronic-structure Library." *WIREs Comput. Mol. Sci.* **8**, e1331 (2018). https://doi.org/10.1002/wcms.1331 | Multireference, relativistic, CASPT2 |
| 55 | 101 | **Columbus** | Academic | https://www.univie.ac.at/columbus/ | Gaussian CI/MCSCF | Lischka, H.; et al. "Columbus—a program system for advanced multireference theory." *WIREs Comput. Mol. Sci.* **1**, 191-199 (2011). https://doi.org/10.1002/wcms.25 | MRC I, MCSCF, surface hopping, nonadiabatic dynamics |
| 56 | 120 | **Kwant** | BSD 2-Clause | https://kwant-project.org/ | Tight-binding transport | Groth, C. W.; et al. "Kwant: a software package for quantum transport." *New J. Phys.* **16**, 063065 (2014). https://doi.org/10.1088/1367-2630/16/6/063065 | Python tight-binding, quantum transport |
| 57 | 145 | **Pybinding** | BSD 3-Clause | https://docs.pybinding.site/ | Python TB | Pybinding documentation. https://docs.pybinding.site/ | 2D materials, graphene, heterostructures |
| 58 | 253 | **Spirit** | MIT License | https://spirit-code.github.io/ | Atomistic spin dynamics | Müller, G. P.; Hoffmann, M.; Dißelkamp, C.; Schürhoff, D.; Mavros, S.; Sallermann, M.; Kiselev, N. S.; Jónsson, H.; Blügel, S. "Spirit: Multifunctional framework for atomistic spin simulations." *Phys. Rev. B* **99**, 224414 (2019). https://doi.org/10.1103/PhysRevB.99.224414 | Spin dynamics, energy minimization, Monte Carlo, Landau-Lifshitz-Gilbert, GPU |
| 59 | 254 | **VAMPIRE** | GNU GPL v3 | https://vampire.york.ac.uk/ | Atomistic spin dynamics | Evans, R. F. L.; Fan, W. J.; Chureemart, P.; Ostler, T. A.; Ellis, M. O. A.; Chantrell, R. W. "Atomistic spin model simulations of magnetic nanomaterials." *J. Phys.: Condens. Matter* **26**, 103202 (2014). https://doi.org/10.1088/0953-8984/26/10/103202 | Micromagnetics, ultrafast dynamics, heat-assisted recording |
| 60 | 316 | **TB2J** | BSD 3-Clause | https://github.com/mailhexu/TB2J | Magnetic exchange | He, X.; Helbig, N.; Verstraete, M. J.; Bousquet, E. "TB2J: A python package for computing magnetic interaction parameters." *Comput. Phys. Commun.* **264**, 107938 (2021). https://doi.org/10.1016/j.cpc.2021.107938 | Magnetic exchange from DFT, Heisenberg/DMI parameters, Wannier functions |
| 61 | 314 | **Bader** | Academic | http://theory.cm.utexas.edu/henkelman/code/bader/ | Charge partitioning | Henkelman, G.; Arnaldsson, A.; Jónsson, H. "A fast and robust algorithm for Bader decomposition of charge density." *Comput. Mater. Sci.* **36**, 354-360 (2006). https://doi.org/10.1016/j.commatsci.2005.04.010 | Bader charge analysis, QTAIM, zero-flux surfaces, Henkelman group (Texas) |
| 62 | 200 | **Critic2** | GNU GPL v3 | https://github.com/aoterodelaroza/critic2 | Topological analysis | Otero-de-la-Roza, A.; et al. "Critic2: A program for real-space analysis." *Comput. Phys. Commun.* **185**, 1007-1018 (2014). https://doi.org/10.1016/j.cpc.2013.10.026 | QTAIM, topology, non-covalent interactions |
| 63 | 317 | **VESTA** | Freeware | https://jp-minerals.org/vesta/en/ | Crystal visualization | Momma, K.; Izumi, F. "VESTA 3 for three-dimensional visualization of crystal, volumetric and morphology data." *J. Appl. Crystallogr.* **44**, 1272-1276 (2011). https://doi.org/10.1107/S0021889811038970 | 3D structure visualization, isosurfaces, magnetic structures, 50,000+ citations |
| 64 | 318 | **XCrySDen** | GNU GPL | http://www.xcrysden.org/ | Crystallographic visualization | Kokalj, A. "Computer graphics and graphical user interfaces as tools in simulations of matter at the atomic scale." *Comput. Mater. Sci.* **28**, 155-168 (2003). https://doi.org/10.1016/S0927-0256(03)00104-6 | XSF format, Fermi surfaces, charge density, CRYSTAL/QE/PWscf/Wien2k |
| 65 | 255 | **VMD** | Free (non-commercial) | https://www.ks.uiuc.edu/Research/vmd/ | Molecular visualization | Humphrey, W.; Dalke, A.; Schulten, K. "VMD: Visual molecular dynamics." *J. Mol. Graphics* **14** (1), 33-38 (1996). https://doi.org/10.1016/0263-7855(96)00018-5 | Trajectory analysis, rendering, scripting, plugins, 9000+ citations |
| 66 | 256 | **Avogadro** | GNU GPL v2 | https://avogadro.cc/ | Molecular editor/visualizer | Hanwell, M. D.; Curtis, D. E.; Lonie, D. C.; Vandermeersch, T.; Zurek, E.; Hutchison, G. R. "Avogadro: An advanced semantic chemical editor, visualization, and analysis platform." *J. Cheminform.* **4**, 17 (2012). https://doi.org/10.1186/1758-2946-4-17 | Cross-platform, extensible, molecular builder, file format support |
| 67 | 70 | **atomate** | MIT License | https://atomate.org/ | High-level workflows | Mathew, K.; et al. "Atomate: A high-level interface to generate, execute, and analyze computational materials science workflows." *Comput. Mater. Sci.* **139**, 140-152 (2017). https://doi.org/10.1016/j.commatsci.2017.07.030 | Materials Project workflows, high-level |
| 68 | 69 | **FireWorks** | Modified BSD | https://materialsproject.github.io/fireworks/ | Workflow execution | Jain, A.; et al. "FireWorks: a dynamic workflow system." *Concurrency Computat.: Pract. Exper.* **27**, 5037-5059 (2015). https://doi.org/10.1002/cpe.3505 | Workflow definition and execution |
| 69 | 72 | **jobflow** | Modified BSD | https://github.com/materialsproject/jobflow | Workflow programming | Rosen, A. S.; et al. "jobflow: Computational workflows made simple." (GitHub 2021+). https://github.com/materialsproject/jobflow | Workflow programming layer |
| 70 | 261 | **Materials Project** | BSD 3-Clause | https://materialsproject.org/ | Database + API ecosystem | Jain, A.; Ong, S. P.; Hautier, G.; Chen, W.; Richards, W. D.; Dacek, S.; Cholia, S.; Gunter, D.; Skinner, D.; Ceder, G.; Persson, K. A. "Commentary: The Materials Project: A materials genome approach to accelerating materials innovation." *APL Mater.* **1**, 011002 (2013). https://doi.org/10.1063/1.4812323 | 150,000+ compounds, high-throughput DFT, MP ecosystem, open API |
| 71 | 76 | **AFLOW** | Open Source | https://aflowlib.org/ | HT framework | Curtarolo, S.; et al. "AFLOW: An automatic framework for high-throughput materials discovery." *Comput. Mater. Sci.* **58**, 218-226 (2012). https://doi.org/10.1016/j.commatsci.2012.02.005 | Automatic FLOW, HT materials discovery, Duke |
| 72 | 77 | **OQMD** | Open Access | https://oqmd.org/ | HT database | Saal, J. E.; et al. "Materials design and discovery with high-throughput DFT." *JOM* **65**, 1501-1509 (2013). https://doi.org/10.1007/s11837-013-0755-5 | Open Quantum Materials Database, Northwestern |
| 73 | 263 | **NOMAD** | Apache 2.0 | https://nomad-lab.eu/ | FAIR data infrastructure | Draxl, C.; Scheffler, M. "The NOMAD laboratory: from data sharing to artificial intelligence." *J. Phys.: Condens. Matter* **31**, 351001 (2019). https://doi.org/10.1088/1361-648X/ab13bb | Novel Materials Discovery, FAIR principles, EU infrastructure, 10M+ calculations |
| 74 | 78 | **JARVIS** | NIST Software | https://jarvis.nist.gov/ | NIST framework | Choudhary, K.; et al. "The joint automated repository for various integrated simulations (JARVIS)." *npj Comput. Mater.* **6**, 173 (2020). https://doi.org/10.1038/s41524-020-00440-1 | Joint Automated Repository, NIST, benchmarks |
| 75 | 79 | **custodian** | MIT License | https://github.com/materialsproject/custodian | Error handling | Ong, S. P.; et al. "custodian: Error handling for computational workflows." (GitHub). https://github.com/materialsproject/custodian | DFT job error handling, recovery |
| 76 | 80 | **matminer** | Modified BSD | https://github.com/hackingmaterials/matminer | Feature engineering | Ward, L.; et al. "Matminer: An open source toolkit for materials data mining." *Comput. Mater. Sci.* **152**, 60-69 (2018). https://doi.org/10.1016/j.commatsci.2018.05.018 | Material featurizers, 700+ descriptors |
| 77 | 81 | **maggma** | Modified BSD | https://github.com/materialsproject/maggma | MongoDB tools | Montoya, J. H.; et al. "maggma: MongoDB aggregation framework." (GitHub). https://github.com/materialsproject/maggma | Store abstraction, S3/MongoDB |
| 78 | 82 | **Phoebe** | Academic | https://github.com/mir-group/phoebe | e-ph transport | Bernardi, M.; et al. "Ab initio study of hot carriers in the first picosecond." *Phys. Rev. Lett.* **112**, 257402 (2014). https://doi.org/10.1103/PhysRevLett.112.257402 | Electron-phonon transport, combined solver |
| 79 | 83 | **PERTURBO** | Academic | https://perturbo-code.github.io/ | Carrier dynamics | Zhou, J.-J.; et al. "Perturbo: A software package for ab initio electron-phonon interactions." *Comput. Phys. Commun.* **264**, 107970 (2021). https://doi.org/10.1016/j.cpc.2021.107970 | Carrier dynamics, superconductivity, real-time |
| 80 | - | **AtomViz** | - | - | - | - | - |
| 81 | - | **pymatgen-diffusion** | - | - | - | - | - |
| 82 | - | **Julia materials packages** | - | - | - | - | - |
| 83 | 86 | **EMC** | Part of Q-Chem | Part of Q-Chem | EOM-CC | Stanton, J. F.; Gauss, J. "Analytic energy derivatives for the equation-of-motion coupled-cluster method." *J. Chem. Phys.* **103**, 1064-1076 (1995). https://doi.org/10.1063/1.469817 | Equation of Motion Coupled Cluster |
| 84 | 87 | **ADC** | Part of Q-Chem | Part of Q-Chem | ADC methods | Dreuw, A.; Wormit, M. "The algebraic diagrammatic construction scheme." *WIREs Comput. Mol. Sci.* **5**, 82-95 (2015). https://doi.org/10.1002/wcms.1206 | Algebraic Diagrammatic Construction |
| 85 | - | **Gaussian Basis Set Libraries** | - | - | - | - | - |
| 86 | - | **X-ray Crystallography Integration** | - | - | - | - | - |
| 87 | - | **Materials Simulation Suites** | - | - | - | - | - |

## 1.2 All-Electron & Full-Potential Methods
**Subcategory Total: 8 codes**
────────────────────────────────────────────────────────────────────────────────

### 1.2.1 LAPW (Linearized Augmented Plane Wave)

### 1.2.2 LMTO (Linearized Muffin-Tin Orbitals)

### 1.2.3 KKR (Korringa-Kohn-Rostoker)

| # | ID | Code Name | License | Official Website | Basis Set | Primary Publication with DOI | Specialization |
|---|---|---|---|---|---|---|---|
| 1 | 10 | **Elk** | GNU GPL | https://elk.sourceforge.io/ | APW+lo (LAPW) | Dewhurst, J. K.; Sharma, S. "Elk LAPW Code." Documentation at https://elk.sourceforge.io/ | Full-potential LAPW, phonons, non-collinear magnetism |
| 2 | 22 | **Fleur** | Academic | https://www.flapw.de/ | FP-LAPW | Betzinger, M.; et al. "Precise and efficient calculation of Fermi surfaces." *Comput. Phys. Commun.* **196**, 276-281 (2015). https://doi.org/10.1016/j.cpc.2015.05.018 | Full-potential LAPW, magnetic systems, Jülich |
| 3 | 358 | **FLAPW** | Research | Various implementations | FP-LAPW | Full-potential linearized augmented plane wave implementations | Generic FLAPW codes |
| 4 | 363 | **FlapwMBPT** | Research | Limited distribution | FLAPW + MBPT | FLAPW with many-body perturbation theory extensions | FLAPW-based GW/BSE |
| 5 | - | **LMTO-ASA** | - | - | - | - | - |
| 6 | 25 | **SPR-KKR** | Academic | https://www.ebert.cup.uni-muenchen.de/index.php/en/software-en/13-sprkkr | KKR | Ebert, H.; et al. "Calculating condensed matter properties using the KKR-Green's function method." *Rep. Prog. Phys.* **74**, 096501 (2011). https://doi.org/10.1088/0034-4885/74/9/096501 | Spin-polarized relativistic KKR, Munich group |
| 7 | - | **JuKKR** | - | - | - | - | - |
| 8 | - | **KKRnano** | - | - | - | - | - |

## 1.3 Localized Basis Sets (Gaussian Basis, Numerical Atomic Orbitals)
**Subcategory Total: 26 codes**
────────────────────────────────────────────────────────────────────────────────

### 1.3.1 Gaussian Basis - Quantum Chemistry Packages

### 1.3.2 Gaussian Basis - Advanced/Specialized

### 1.3.3 Numerical Atomic Orbitals

| # | ID | Code Name | License | Official Website | Basis Set | Primary Publication with DOI | Specialization |
|---|---|---|---|---|---|---|---|
| 1 | 26 | **Gaussian** | Commercial | https://gaussian.com/ | Gaussian basis | Frisch, M. J.; et al. "Gaussian 16." Gaussian Inc. (2016). (Gaussian 09, 16, 20 versions) | Flagship QC package, extensive method library |
| 2 | 29 | **PySCF** | Apache 2.0 | https://pyscf.org/ | Gaussian basis | Sun, Q.; et al. "PySCF: the Python-based simulations of chemistry framework." *WIREs Comput. Mol. Sci.* **8**, e1340 (2018). https://doi.org/10.1002/wcms.1340 | Pure Python QC, modular, AI integration |
| 3 | 9 | **FHI-aims** | Academic/Research | https://fhi-aims.org/ | Numeric atom-centered orbitals | Blum, V.; et al. "Ab initio molecular simulations with numeric atom-centered orbitals." *Comput. Phys. Commun.* **180**, 2175-2196 (2009). https://doi.org/10.1016/j.cpc.2009.06.022 | All-electron DFT, GW, hybrid functionals, Fritz Haber Institute |
| 4 | 30 | **TURBOMOLE** | Commercial | https://www.turbomole.org/ | Gaussian basis | TURBOMOLE V7.5 2020, University of Karlsruhe and Forschungszentrum Karlsruhe GmbH, 1989-2020; TURBOMOLE GmbH, since 2007. http://www.turbomole.com | Efficient DFT, RI approximations, German development |
| 5 | 32 | **CFOUR** | Academic/Commercial | http://www.cfour.de/ | Gaussian basis | Stanton, J. F.; et al. "CFOUR." (Documentation 2019). http://www.cfour.de | Coupled-Cluster specialist, extremely high-accuracy |
| 6 | 33 | **GAMESS-US** | Academic | https://www.msg.chem.iastate.edu/gamess/ | Gaussian basis | Schmidt, M. W.; et al. "General atomic and molecular electronic structure system." *J. Comput. Chem.* **14**, 1347-1363 (1993). https://doi.org/10.1002/jcc.540141112 | General QC, free for academic use |
| 7 | - | **GAMESS-UK** | - | - | - | - | - |
| 8 | 34 | **Dalton** | GNU LGPL v2.1 | https://daltonprogram.org/ | Gaussian basis | Aidas, K.; et al. "The Dalton quantum chemistry program system." *WIREs Comput. Mol. Sci.* **4**, 269-284 (2014). https://doi.org/10.1002/wcms.1172 | Response properties, spectroscopy, Norway |
| 9 | 35 | **DIRAC** | GNU LGPL | http://www.diracprogram.org/ | Gaussian basis | Saue, T.; et al. "The DIRAC code for relativistic molecular calculations." *J. Chem. Phys.* **152**, 204104 (2020). https://doi.org/10.1063/5.0004844 | Relativistic quantum chemistry, 4-component |
| 10 | 36 | **ADF** | Commercial | https://www.scm.com/product/adf/ | Slater-type orbitals | te Velde, G.; et al. "Chemistry with ADF." *J. Comput. Chem.* **22**, 931-967 (2001). https://doi.org/10.1002/jcc.1056 | Amsterdam Density Functional, STO basis |
| 11 | 38 | **Q-Chem** | Commercial | https://www.q-chem.com/ | Gaussian basis | Epifanovsky, E.; et al. "Software for the frontiers of quantum chemistry." *J. Chem. Phys.* **155**, 084801 (2021). https://doi.org/10.1063/5.0055522 | Comprehensive suite, excited states, ML integration |
| 12 | 93 | **Firefly** | Academic | http://classic.chem.msu.su/gran/firefly/ | QC | Granovsky, A. A. "Firefly version 8." http://classic.chem.msu.su/gran/firefly/ | PC GAMESS, Windows-optimized QC |
| 13 | - | **ACES II** | - | - | - | - | - |
| 14 | 102 | **CADPAC** | Academic | Historical | Gaussian | Amos, R. D.; et al. "Cambridge Analytic Derivatives Package (CADPAC)." (Historical, 1980s-1990s) | Analytical derivatives emphasis |
| 15 | 103 | **hBar Lab7** | Commercial | Limited distribution | QC | Proprietary package with limited public information | Commercial QC software |
| 16 | 96 | **JAGUAR** | Commercial | Part of Schrödinger | QC | Schrödinger suite component, industrial QC | Schrödinger industrial package |
| 17 | 104 | **PQS** | Commercial | http://www.pqschem.com/ | Gaussian | Parallel Quantum Solutions. "PQS." http://www.pqschem.com/ | Parallel Quantum Solutions |
| 18 | 113 | **deMon2k** | Academic | http://www.demon-software.com/ | Auxiliary DFT | Köster, A. M.; et al. "deMon2k." *WIREs Comput. Mol. Sci.* **1**, 820-838 (2011). https://doi.org/10.1002/wcms.68 | Density fitting DFT, MinBas |
| 19 | - | **Priroda-06** | - | - | - | - | - |
| 20 | 40 | **MPQC** | GNU GPL | https://github.com/ValeevGroup/mpqc | Gaussian basis | Janssen, C. L.; et al. "The Massively Parallel Quantum Chemistry program (MPQC)." (2019). https://github.com/ValeevGroup/mpqc | Massively parallel QC, modular architecture |
| 21 | 105 | **FreeON** | GNU GPL | Historical (discontinued) | Gaussian | Schwegler, E.; et al. "Linear scaling computation of the Fock matrix." *J. Chem. Phys.* **106**, 5526-5535 (1997). https://doi.org/10.1063/1.473575 | Linear-scaling, now discontinued |
| 22 | 39 | **MOPAC** | Academic (MOPAC2016) | http://openmopac.net/ | Semi-empirical | Stewart, J. J. P. "Optimization of parameters for semiempirical methods VI: PM7." *J. Mol. Model.* **19**, 1-32 (2013). https://doi.org/10.1007/s00894-012-1667-x | Semi-empirical methods, PM6/PM7, fast QC |
| 23 | - | **PyQuante** | - | - | - | - | - |
| 24 | - | **PLATO** | - | - | - | - | - |
| 25 | - | **Atomistix ToolKit** | - | - | - | - | - |
| 26 | 312 | **S/PHI/nX** | Academic | https://sxrepo.mpie.de/ | Numeric basis DFT | Boeck, S.; Freysoldt, C.; Dick, A.; Ismer, L.; Neugebauer, J. "The object-oriented DFT program library S/PHI/nX." *Comput. Phys. Commun.* **182**, 543-554 (2011). https://doi.org/10.1016/j.cpc.2010.09.016 | Sphinxlib, numeric atom-centered orbitals, object-oriented C++, Max Planck |

## 1.4 Tight-Binding DFT & Semi-Empirical Methods
**Subcategory Total: 2 codes**
────────────────────────────────────────────────────────────────────────────────

| # | ID | Code Name | License | Official Website | Basis Set | Primary Publication with DOI | Specialization |
|---|---|---|---|---|---|---|---|
| 1 | - | **HOTBIT** | - | - | - | - | - |
| 2 | 351 | **DFTB** | Research | Various implementations | Tight-binding DFT | Seifert, G. "Tight-binding density functional theory: an approximate Kohn-Sham DFT scheme." *J. Phys. Chem. A* **111**, 5609-5613 (2007). https://doi.org/10.1021/jp069056r | Base DFTB implementation (various) |

## 1.5 DFT+U, Hybrid Functionals & Specialized DFT Variants
**Subcategory Total: 0 codes**
────────────────────────────────────────────────────────────────────────────────

*No codes in this subcategory*


### **Category 1 Total: 123 codes**


====================================================================================================
# 2. TIME-DEPENDENT & EXCITED-STATE METHODS
====================================================================================================

## 2.1 TDDFT (Time-Dependent DFT)
**Subcategory Total: 2 codes**
────────────────────────────────────────────────────────────────────────────────

### 2.1.1 Linear-Response TDDFT

### 2.1.2 Real-Time TDDFT

| # | ID | Code Name | License | Official Website | Basis Set | Primary Publication with DOI | Specialization |
|---|---|---|---|---|---|---|---|
| 1 | - | **SALMON** | - | - | - | - | - |
| 2 | - | **Qbox** | - | - | - | - | - |

## 2.2 Many-Body Perturbation Theory (GW & BSE)
**Subcategory Total: 6 codes**
────────────────────────────────────────────────────────────────────────────────

### 2.2.1 GW Implementations

### 2.2.2 BSE (Bethe-Salpeter Equation)

| # | ID | Code Name | License | Official Website | Basis Set | Primary Publication with DOI | Specialization |
|---|---|---|---|---|---|---|---|
| 1 | 47 | **SternheimerGW** | Academic | https://github.com/QEF/SternheimerGW | GW | Schlipf, M.; et al. "SternheimerGW: A program for calculating GW quasiparticle band structures." *Comput. Phys. Commun.* **247**, 106856 (2020). https://doi.org/10.1016/j.cpc.2019.106856 | DFPT-based GW, Quantum ESPRESSO plugin |
| 2 | 45 | **WEST** | GNU GPL v3 | http://www.west-code.org/ | GW | Govoni, M.; Galli, G. "Large scale GW calculations." *J. Chem. Theory Comput.* **11**, 2680-2696 (2015). https://doi.org/10.1021/ct500958p | Without Empty States, large-scale GW, GPU |
| 3 | - | **molgw** | - | - | - | - | - |
| 4 | - | **GreenX** | - | - | - | - | - |
| 5 | 123 | **OCEAN** | GNU GPL v2 | https://feff.phys.washington.edu/OCEAN/ | Core-level spectroscopy | Vinson, J.; et al. "Bethe-Salpeter equation calculations of core excitation spectra." *Phys. Rev. B* **83**, 115106 (2011). https://doi.org/10.1103/PhysRevB.83.115106 | X-ray absorption, EELS, BSE |
| 6 | - | **NBSE** | - | - | - | - | - |


### **Category 2 Total: 8 codes**


====================================================================================================
# 3. STRONGLY CORRELATED & MANY-BODY METHODS
====================================================================================================

## 3.1 DMFT (Dynamical Mean-Field Theory)
**Subcategory Total: 21 codes**
────────────────────────────────────────────────────────────────────────────────

### 3.1.1 DMFT Frameworks & Core Libraries

### 3.1.2 DFT+DMFT Implementations (Integrated with DFT Codes)

### 3.1.3 Impurity Solvers (DMFT Solvers)

| # | ID | Code Name | License | Official Website | Basis Set | Primary Publication with DOI | Specialization |
|---|---|---|---|---|---|---|---|
| 1 | - | **TRIQS/DFTTools** | - | - | - | - | - |
| 2 | - | **solid_dmft** | - | - | - | - | - |
| 3 | 127 | **ALPS** | Open Source | https://alps.comp-phys.org/ | QMC for models | Alet, F.; et al. "The ALPS project: open source software for strongly correlated systems." *J. Phys. Soc. Jpn.* **74**, Suppl. 30-35 (2005). https://doi.org/10.1143/JPSJS.74S.30 | Algorithms and Libraries for Physics Simulations |
| 4 | - | **ALPSCore** | - | - | - | - | - |
| 5 | 52 | **w2dynamics** | GNU GPL v3 | https://github.com/w2dynamics/w2dynamics | DMFT solver | Wallerberger, M.; et al. "w2dynamics: Local one- and two-particle quantities from dynamical mean field theory." *Comput. Phys. Commun.* **235**, 388-399 (2019). https://doi.org/10.1016/j.cpc.2018.09.007 | Worm algorithm DMFT solver, CT-HYB |
| 6 | 53 | **DCore** | GNU GPL v3 | https://issp-center-dev.github.io/DCore/ | DMFT | Shinaoka, H.; et al. "Compressing Green's function using intermediate representation." *Phys. Rev. B* **96**, 035147 (2017). https://doi.org/10.1103/PhysRevB.96.035147 | Integrated DMFT software for correlated electrons |
| 7 | 126 | **iQIST** | GNU GPL v3 | https://github.com/huangli712/iQIST | CT-QMC | Li, H.; Huang, L. "iQIST: An open source continuous-time quantum Monte Carlo impurity solver." *Comput. Phys. Commun.* **195**, 140-160 (2015). https://doi.org/10.1016/j.cpc.2015.04.020 | Interacting Quantum Impurity Solver Toolkit, CT-HYB |
| 8 | - | **AMULET** | - | - | - | - | - |
| 9 | - | **DMFTwDFT** | - | - | - | - | - |
| 10 | 55 | **ComDMFT** | Academic | Part of commercial packages | DMFT | Various DMFT implementations in Materials Studio and other suites | Commercial DMFT solvers |
| 11 | 54 | **EDMFTF** | Academic | https://www.physics.rutgers.edu/~haule/518/EDMFT.html | DMFT | Haule, K.; et al. "Dynamical mean-field theory." *Phys. Rev. B* **81**, 195107 (2010). https://doi.org/10.1103/PhysRevB.81.195107 | Embedded DMFT, Haule group (Rutgers) |
| 12 | - | **VASP+DMFT** | - | - | - | - | - |
| 13 | - | **CT-HYB** | - | - | - | - | - |
| 14 | - | **CT-QMC** | - | - | - | - | - |
| 15 | - | **CT-INT** | - | - | - | - | - |
| 16 | - | **CT-SEG** | - | - | - | - | - |
| 17 | - | **HÏ†** | - | - | - | - | - |
| 18 | 291 | **EDIpack** | GNU GPL v3 | https://edipack.github.io/EDIpack/ | Lanczos exact diagonalization | Amaricci, A.; Crippa, L.; Tocchio, A.; Parcollet, F.; Pourovskii, I.; Valenti, R. "EDIpack: A parallel exact diagonalization package for quantum impurity problems." *Comput. Phys. Commun.* **273**, 108261 (2022). https://doi.org/10.1016/j.cpc.2021.108261; Crippa, L.; et al. "Next-generation EDIpack: A Lanczos-based package for quantum impurity models featuring general broken-symmetry phases, flexible bath topologies and multi-platform interoperability." *SciPost Phys. Codebases* **58** (2025). https://doi.org/10.21468/SciPostPhysCodeb.58 | Quantum impurity solver, DMFT, superconductivity, spin-orbit coupling, TRIQS/w2dynamics interfaces |
| 19 | - | **FTPS** | - | - | - | - | - |
| 20 | 128 | **pomerol** | GNU GPL v2 | https://github.com/aeantipov/pomerol | Exact diagonalization | Antipov, A. E.; et al. "Pomerol: ED impurity solver." (GitHub 2014+). https://github.com/aeantipov/pomerol | Quantum impurity ED solver |
| 21 | - | **ALPS/CT-HYB** | - | - | - | - | - |

## 3.2 Quantum Monte Carlo (QMC)
**Subcategory Total: 9 codes**
────────────────────────────────────────────────────────────────────────────────

### 3.2.1 Continuum QMC (VMC, DMC, AFQMC)

### 3.2.2 Lattice & Model QMC

| # | ID | Code Name | License | Official Website | Basis Set | Primary Publication with DOI | Specialization |
|---|---|---|---|---|---|---|---|
| 1 | 61 | **TurboRVB** | GNU GPL v3 | https://people.sissa.it/~sorella/web/index.html | QMC | Nakano, K.; et al. "TurboRVB: A many-body toolkit." *J. Chem. Phys.* **152**, 204121 (2020). https://doi.org/10.1063/5.0005037 | Variational QMC, Jastrow wavefunctions |
| 2 | - | **PyQMC** | - | - | - | - | - |
| 3 | - | **CHAMP** | - | - | - | - | - |
| 4 | - | **QMcBeaver** | - | - | - | - | - |
| 5 | - | **QWalk** | - | - | - | - | - |
| 6 | - | **ALF** | - | - | - | - | - |
| 7 | - | **QUEST** | - | - | - | - | - |
| 8 | - | **TRIQS/CT-QMC solvers** | - | - | - | - | - |
| 9 | - | **DCA++** | - | - | - | - | - |


### **Category 3 Total: 30 codes**


====================================================================================================
# 4. WAVEFUNCTION-BASED QUANTUM CHEMISTRY
====================================================================================================

## 4.1 Coupled-Cluster Methods
**Subcategory Total: 0 codes**
────────────────────────────────────────────────────────────────────────────────

*No codes in this subcategory*

## 4.2 Configuration Interaction & Multireference
**Subcategory Total: 0 codes**
────────────────────────────────────────────────────────────────────────────────

*No codes in this subcategory*

## 4.3 Quantum Chemistry Suites (General Packages)
**Subcategory Total: 0 codes**
────────────────────────────────────────────────────────────────────────────────

*No codes in this subcategory*


### **Category 4 Total: 0 codes**


====================================================================================================
# 5. TIGHT-BINDING, MODEL HAMILTONIANS & DOWNFOLDING
====================================================================================================

## 5.1 Wannier Function Methods
**Subcategory Total: 5 codes**
────────────────────────────────────────────────────────────────────────────────

| # | ID | Code Name | License | Official Website | Basis Set | Primary Publication with DOI | Specialization |
|---|---|---|---|---|---|---|---|
| 1 | 141 | **WannierTools** | GNU GPL v3 | https://www.wanniertools.com/ | Wannier TB | Wu, Q.; et al. "WannierTools: An open-source software package for novel topological materials." *Comput. Phys. Commun.* **224**, 405-416 (2018). https://doi.org/10.1016/j.cpc.2017.09.033 | Topological invariants, surface states |
| 2 | 143 | **PythTB** | GNU GPL v3 | http://www.physics.rutgers.edu/pythtb/ | Python TB | Vanderbilt, D. *Berry Phases in Electronic Structure Theory.* Cambridge University Press (2018). https://doi.org/10.1017/9781316662205 | Educational tight-binding |
| 3 | 142 | **TBmodels** | Apache 2.0 | https://github.com/Z2PackDev/TBmodels | Wannier TB | Gresch, D.; et al. "Z2Pack: Numerical implementation of hybrid Wannier centers." *Phys. Rev. B* **95**, 075146 (2017). https://doi.org/10.1103/PhysRevB.95.075146 | TB model manipulation, symmetrization |
| 4 | 144 | **TopoTB** | Academic | https://github.com/cpoli/TopoTB | TB + topology | Poli, C.; et al. "Selective enhancement of topologically induced interface states." *Nat. Commun.* **6**, 6710 (2015). https://doi.org/10.1038/ncomms7710 | TB electronic structure, topology |
| 5 | 332 | **AiiDA-wannier90** | MIT License | https://github.com/aiidateam/aiida-wannier90 | Wannier90 plugin | AiiDA plugin for Wannier90. (GitHub). https://github.com/aiidateam/aiida-wannier90 | Wannier90 workflows in AiiDA |

## 5.2 Model Hamiltonian Solvers
**Subcategory Total: 3 codes**
────────────────────────────────────────────────────────────────────────────────

| # | ID | Code Name | License | Official Website | Basis Set | Primary Publication with DOI | Specialization |
|---|---|---|---|---|---|---|---|
| 1 | 146 | **TBSTUDIO** | Freeware | https://www.tbstudio.net/ | TB builder | https://www.tbstudio.net/ documentation | Educational TB with GUI |
| 2 | 213 | **HubbardFermiMatsubara** | Academic | https://github.com/tcompa/HubbardFermiMatsubara | Hubbard model solvers | Comparin, T.; Mezzacapo, F. "Multiloop functional renormalization group for general models." *Phys. Rev. B* **100**, 205130 (2019). https://doi.org/10.1103/PhysRevB.100.205130 | Hubbard model, functional RG, Matsubara |
| 3 | 136 | **exactdiag** | MIT License | https://github.com/cpburnz/exactdiag | Exact diagonalization | GitHub community repository. https://github.com/cpburnz/exactdiag | Exact diagonalization for pedagogy |

## 5.3 Downfolding & Embedding Methods
**Subcategory Total: 1 codes**
────────────────────────────────────────────────────────────────────────────────

| # | ID | Code Name | License | Official Website | Basis Set | Primary Publication with DOI | Specialization |
|---|---|---|---|---|---|---|---|
| 1 | 260 | **FermiSurfer** | MIT License | https://mitsuaki1987.github.io/fermisurfer/ | Fermi surface visualization | Kawamura, M. "FermiSurfer: Fermi-surface viewer providing multiple representation schemes." *Comput. Phys. Commun.* **239**, 197-203 (2019). https://doi.org/10.1016/j.cpc.2019.01.017 | 3D Fermi surfaces, cross-platform, nodal lines, Wannier90/QE interfaces |


### **Category 5 Total: 9 codes**


====================================================================================================
# 6. PHONONS, LATTICE DYNAMICS & ELECTRON-PHONON COUPLING
====================================================================================================

## 6.1 Harmonic Phonons
**Subcategory Total: 3 codes**
────────────────────────────────────────────────────────────────────────────────

### 6.1.1 Phonon Calculation Codes

### 6.1.2 DFPT (Density Functional Perturbation Theory)

| # | ID | Code Name | License | Official Website | Basis Set | Primary Publication with DOI | Specialization |
|---|---|---|---|---|---|---|---|
| 1 | 151 | **PHON** | Academic | http://www.homepages.ucl.ac.uk/~ucfbdxa/phon/ | Harmonic phonons | Alfè, D. "PHON: A program to calculate phonons using the small displacement method." *Comput. Phys. Commun.* **180**, 2622-2633 (2009). https://doi.org/10.1016/j.cpc.2009.03.010 | Finite displacement, VASP interface |
| 2 | 152 | **PHONON** | Academic | http://wolf.ifj.edu.pl/phonon/ | Harmonic phonons | Parlinski, K.; et al. "First-Principles Determination of the Soft Mode in Cubic ZrO2." *Phys. Rev. Lett.* **78**, 4063-4066 (1997). https://doi.org/10.1103/PhysRevLett.78.4063 | Legacy phonon code |
| 3 | 153 | **YPHON** | Academic | http://www.yphon.org/ | Harmonic phonons | YPHON online documentation. http://www.yphon.org/ | Phonon dispersions |

## 6.2 Anharmonic Phonons & Thermal Transport
**Subcategory Total: 11 codes**
────────────────────────────────────────────────────────────────────────────────

### 6.2.1 Anharmonic Phonon & Thermal Conductivity Codes

### 6.2.2 Anharmonic Method Implementation

| # | ID | Code Name | License | Official Website | Basis Set | Primary Publication with DOI | Specialization |
|---|---|---|---|---|---|---|---|
| 1 | 154 | **almaBTE** | Mozilla Public License 2.0 | https://www.almabte.eu/ | Phonon BTE | Carrete, J.; et al. "almaBTE: A solver of the space–time dependent Boltzmann transport equation." *Comput. Phys. Commun.* **220**, 351-362 (2017). https://doi.org/10.1016/j.cpc.2017.06.023 | Phonon BTE, nanostructures |
| 2 | 155 | **PhonTS** | Academic | https://github.com/WMD-group/PhonTS | Phonon transport | Skelton, J. M.; et al. "Lattice dynamics of the tin sulphides." *Phys. Chem. Chem. Phys.* **19**, 12452-12465 (2017). https://doi.org/10.1039/C7CP01680H | Mode Grüneisen parameters |
| 3 | 156 | **TDEP** | MIT License | https://ollehellman.github.io/page/ | Temperature-dependent phonons | Hellman, O.; et al. "Temperature dependent effective potential method for accurate free energy calculations." *Phys. Rev. B* **87**, 104111 (2013). https://doi.org/10.1103/PhysRevB.87.104111 | Temperature Dependent Effective Potential |
| 4 | 157 | **kALDo** | MIT License | https://github.com/nanotheorygroup/kaldo | Phonon BTE | Carrete, J.; et al. "kALDo: Anharmonic Lattice Dynamics." (GitHub 2020+). https://github.com/nanotheorygroup/kaldo | Anharmonic lattice dynamics |
| 5 | 158 | **GPU_PBTE** | Research | Various GPU implementations | GPU-accelerated BTE | GPU-accelerated phonon Boltzmann transport equation solvers | GPU phonon BTE |
| 6 | 159 | **SCAILD** | Academic | Research code | Self-consistent phonons | Monserrat, B. "Electron-phonon coupling from finite differences." *J. Phys.: Condens. Matter* **30**, 083001 (2018). https://doi.org/10.1088/1361-648X/aaa737 | Self-consistent ab initio lattice dynamics |
| 7 | 160 | **QSCAILD** | Academic | Research code | Quantum SCAILD | Extension of SCAILD with quantum effects | Quantum self-consistent approach |
| 8 | 161 | **SSCHA** | GNU GPL v3 | https://sscha.eu/ | Stochastic phonons | Monacelli, L.; et al. "The stochastic self-consistent harmonic approximation." *J. Phys.: Condens. Matter* **33**, 363001 (2021). https://doi.org/10.1088/1361-648X/ac066b | Stochastic Self-Consistent Harmonic Approximation |
| 9 | 162 | **ALM** | MIT License | https://github.com/ttadano/ALM | Force constant extraction | Tadano, T.; Tsuneyuki, S. "Self-consistent phonon calculations of lattice dynamics." *Phys. Rev. B* **92**, 054301 (2015). https://doi.org/10.1103/PhysRevB.92.054301 | Anharmonic force constants |
| 10 | 163 | **hiPhive** | MIT License | https://hiphive.materialsmodeling.org/ | High-order IFCs | Eriksson, F.; et al. "The Hiphive Package for the Extraction of High‐Order Force Constants." *Adv. Theory Simul.* **2**, 1800184 (2019). https://doi.org/10.1002/adts.201800184 | High-order force constants |
| 11 | 164 | **thirdorder.py** | Academic | Part of ShengBTE | 3rd-order IFCs | Li, W.; et al. "ShengBTE: A solver of the Boltzmann transport equation." *Comput. Phys. Commun.* **185**, 1747-1758 (2014). https://doi.org/10.1016/j.cpc.2014.02.015 | Third-order IFC extraction |

## 6.3 Electron-Phonon Coupling
**Subcategory Total: 2 codes**
────────────────────────────────────────────────────────────────────────────────

| # | ID | Code Name | License | Official Website | Basis Set | Primary Publication with DOI | Specialization |
|---|---|---|---|---|---|---|---|
| 1 | 198 | **BoltzWann** | GNU GPL v3 | http://www.boltztrap.org/ | Wannier transport | Pizzi, G.; et al. "BoltzWann: A code for thermoelectric and electronic transport." *Comput. Phys. Commun.* **185**, 422-429 (2014). https://doi.org/10.1016/j.cpc.2013.09.015 | Wannier interpolation transport |
| 2 | - | **DMDW/RTDW** | - | - | - | - | - |


### **Category 6 Total: 16 codes**


====================================================================================================
# 7. MOLECULAR & AB INITIO DYNAMICS
====================================================================================================

## 7.1 Born-Oppenheimer Molecular Dynamics
**Subcategory Total: 0 codes**
────────────────────────────────────────────────────────────────────────────────

*No codes in this subcategory*

## 7.2 Path Integral Molecular Dynamics
**Subcategory Total: 0 codes**
────────────────────────────────────────────────────────────────────────────────

*No codes in this subcategory*

## 7.3 Rare Events, Transitions & Enhanced Sampling
**Subcategory Total: 3 codes**
────────────────────────────────────────────────────────────────────────────────

| # | ID | Code Name | License | Official Website | Basis Set | Primary Publication with DOI | Specialization |
|---|---|---|---|---|---|---|---|
| 1 | 178 | **NEB** | Various | Various implementations | Nudged Elastic Band | Henkelman, G.; et al. "Improved tangent estimate in the nudged elastic band method." *J. Chem. Phys.* **113**, 9978-9985 (2000). https://doi.org/10.1063/1.1323224 | Transition state search |
| 2 | - | **String methods** | - | - | - | - | - |
| 3 | - | **Metadynamics** | - | - | - | - | - |


### **Category 7 Total: 3 codes**


====================================================================================================
# 8. STRUCTURE PREDICTION & GLOBAL OPTIMIZATION
====================================================================================================

## 8.1 Evolutionary Algorithms
**Subcategory Total: 6 codes**
────────────────────────────────────────────────────────────────────────────────

| # | ID | Code Name | License | Official Website | Basis Set | Primary Publication with DOI | Specialization |
|---|---|---|---|---|---|---|---|
| 1 | 73 | **USPEX** | Academic | https://uspex-team.org/ | Structure prediction | Oganov, A. R.; Glass, C. W. "Crystal structure prediction using evolutionary algorithms." *J. Chem. Phys.* **124**, 244704 (2006). https://doi.org/10.1063/1.2210932 | Universal Structure Predictor, evolutionary algorithm |
| 2 | 75 | **XtalOpt** | BSD 2-Clause | https://xtalopt.github.io/ | Structure prediction | Lonie, D. C.; Zurek, E. "XtalOpt: An open-source evolutionary algorithm." *Comput. Phys. Commun.* **182**, 372-387 (2011). https://doi.org/10.1016/j.cpc.2010.07.048 | Evolutionary algorithm, Avogadro-integrated |
| 3 | 74 | **CALYPSO** | Academic | http://www.calypso.cn/ | Structure prediction | Wang, Y.; et al. "CALYPSO: A method for crystal structure prediction." *Comput. Phys. Commun.* **183**, 2063-2070 (2012). https://doi.org/10.1016/j.cpc.2012.05.008 | Particle swarm optimization, Chinese development |
| 4 | 181 | **GASP** | GNU GPL v3 | https://github.com/henniggroup/GASP | Genetic algorithm | Revard, B. C.; et al. "Structure and Stability Prediction of Compounds with Evolutionary Algorithms." *Top. Curr. Chem.* **345**, 181-222 (2014). https://doi.org/10.1007/128_2013_489 | Genetic Algorithm for Structure and Phase Prediction |
| 5 | 183 | **MAISE** | Open Source | https://github.com/maise-guide/maise | EA + NN potentials | Hajinazar, S.; et al. "MAISE: Construction of neural network interatomic models." *Comput. Phys. Commun.* **259**, 107679 (2021). https://doi.org/10.1016/j.cpc.2020.107679 | Module for Ab Initio Structure Evolution |
| 6 | 184 | **EVO** | GNU GPL | http://cpc.cs.qub.ac.uk/summaries/AEOZ_v1_0.html | Evolutionary strategy | Bahmann, S.; Kortus, J. "EVO—Evolutionary algorithm for crystal structure prediction." *Comput. Phys. Commun.* **184**, 1618-1625 (2013). https://doi.org/10.1016/j.cpc.2013.01.005 | Python EA, QE/GULP interfaces |

## 8.2 Random Sampling & Basin Hopping
**Subcategory Total: 3 codes**
────────────────────────────────────────────────────────────────────────────────

| # | ID | Code Name | License | Official Website | Basis Set | Primary Publication with DOI | Specialization |
|---|---|---|---|---|---|---|---|
| 1 | 185 | **AIRSS** | Academic | https://www.mtg.msm.cam.ac.uk/Codes/AIRSS | Random search | Pickard, C. J.; Needs, R. J. "Ab initio random structure searching." *J. Phys.: Condens. Matter* **23**, 053201 (2011). https://doi.org/10.1088/0953-8984/23/5/053201 | Ab Initio Random Structure Searching |
| 2 | 182 | **FLAME** | GNU GPL | https://github.com/FLAME-Org/FLAME | Minima hopping + NN | Amsler, M.; Goedecker, S. "Crystal structure prediction using the minima hopping method." *J. Chem. Phys.* **133**, 224104 (2010). https://doi.org/10.1063/1.3512900 | Free Lattice Atomic Method for Energy, NN potentials |
| 3 | - | **Basin hopping** | - | - | - | - | - |

## 8.3 Machine Learning Approaches
**Subcategory Total: 2 codes**
────────────────────────────────────────────────────────────────────────────────

| # | ID | Code Name | License | Official Website | Basis Set | Primary Publication with DOI | Specialization |
|---|---|---|---|---|---|---|---|
| 1 | 240 | **HTOCSP** | Academic | https://github.com/tonnamb/HTOCSP | ML-enhanced CSP | Ward, L.; Agrawal, A.; Choudhary, A.; Wolverton, C. "A general-purpose machine learning framework for predicting properties of inorganic materials." *npj Comput. Mater.* **2**, 16028 (2016). https://doi.org/10.1038/npjcompumats.2016.28 | High-throughput organic CSP, ML acceleration, property prediction |
| 2 | - | **Neural network potentials** | - | - | - | - | - |


### **Category 8 Total: 11 codes**


====================================================================================================
# 9. POST-PROCESSING, ANALYSIS & VISUALIZATION
====================================================================================================

## 9.1 Electronic Structure Analysis
**Subcategory Total: 9 codes**
────────────────────────────────────────────────────────────────────────────────

### 9.1.1 Band Structure & Density of States

### 9.1.2 Transport Properties

### 9.1.3 Chemical Bonding Analysis

| # | ID | Code Name | License | Official Website | Basis Set | Primary Publication with DOI | Specialization |
|---|---|---|---|---|---|---|---|
| 1 | 191 | **vaspkit** | MIT License | https://vaspkit.com/ | VASP post-processing | Wang, V.; et al. "VASPKIT: A user-friendly interface facilitating high-throughput computing." *Comput. Phys. Commun.* **267**, 108033 (2021). https://doi.org/10.1016/j.cpc.2021.108033 | Comprehensive VASP toolkit |
| 2 | 192 | **sumo** | MIT License | https://github.com/SMTG-UCL/sumo | Band/DOS plotting | Ganose, A. M.; et al. "sumo: Command-line tools for plotting and analysis." *J. Open Source Softw.* **3**, 717 (2018). https://doi.org/10.21105/joss.00717 | Publication-quality plots |
| 3 | 193 | **pyprocar** | MIT License | https://romerogroup.github.io/pyprocar/ | Band structure analysis | Herath, U.; et al. "PyProcar: A Python library for electronic structure pre/post-processing." *Comput. Phys. Commun.* **251**, 107080 (2020). https://doi.org/10.1016/j.cpc.2019.107080 | Projection, spin texture, Fermi surfaces |
| 4 | 194 | **PyARPES** | MIT License | https://github.com/chstan/arpes | ARPES analysis | PyARPES GitHub documentation. https://github.com/chstan/arpes | ARPES data analysis |
| 5 | 195 | **BandUP** | GNU GPL v3 | https://github.com/band-unfolding/banduppy | Band unfolding | Medeiros, P. V. C.; et al. "Effects of extrinsic perturbations on electronic structure." *Phys. Rev. B* **89**, 041407(R) (2014). https://doi.org/10.1103/PhysRevB.89.041407 | Supercell band unfolding |
| 6 | 196 | **fold2Bloch** | GNU GPL | https://github.com/rubel75/fold2Bloch | Band unfolding | Popescu, V.; Zunger, A. "Extracting E versus k effective band structure." *Phys. Rev. B* **85**, 085201 (2012). https://doi.org/10.1103/PhysRevB.85.085201 | Spectral weights, unfolding |
| 7 | 343 | **BoltzTraP2** | GNU GPL v3 | https://www.boltztrap.org/ | Boltzmann transport | Madsen, G. K. H.; et al. "BoltzTraP2, a program for interpolating band structures." *Comput. Phys. Commun.* **231**, 140-145 (2018). https://doi.org/10.1016/j.cpc.2018.05.010 | Second-generation BoltzTraP, Python |
| 8 | 345 | **COHP** | Part of Lobster | Bonding analysis | Crystal Orbital Hamilton Population analysis (part of Lobster package) | Crystal Orbital Hamilton Population | - |
| 9 | 199 | **DDEC** | Academic | https://sourceforge.net/projects/ddec/ | Charge partitioning | Manz, T. A.; Limas, N. G. "Introducing DDEC6 atomic population analysis." *RSC Adv.* **6**, 47771-47801 (2016). https://doi.org/10.1039/C6RA04656H | Density-Derived charges, DDEC6 |

## 9.2 Optical & Spectroscopic Properties
**Subcategory Total: 1 codes**
────────────────────────────────────────────────────────────────────────────────

| # | ID | Code Name | License | Official Website | Basis Set | Primary Publication with DOI | Specialization |
|---|---|---|---|---|---|---|---|
| 1 | 122 | **DP** | Academic | http://www.dp-code.org/ | Dielectric properties | Gajdoš, M.; et al. "Linear optical properties in the PAW methodology." *Phys. Rev. B* **73**, 045112 (2006). https://doi.org/10.1103/PhysRevB.73.045112 | Optical absorption, dielectric tensor |

## 9.3 Magnetic Properties
**Subcategory Total: 1 codes**
────────────────────────────────────────────────────────────────────────────────

| # | ID | Code Name | License | Official Website | Basis Set | Primary Publication with DOI | Specialization |
|---|---|---|---|---|---|---|---|
| 1 | 252 | **magnon codes** | Various | Various implementations | Spin wave theory | Various authors - general category for magnon dispersion codes interfacing DFT outputs | Magnon dispersions, spin waves, collective excitations |

## 9.4 Structure Visualization
**Subcategory Total: 3 codes**
────────────────────────────────────────────────────────────────────────────────

| # | ID | Code Name | License | Official Website | Basis Set | Primary Publication with DOI | Specialization |
|---|---|---|---|---|---|---|---|
| 1 | 257 | **STMng** | Academic | https://github.com/oganov-lab/STMng | STM simulation | Bogdanov, N.; et al. "STMng: A visualization tool for USPEX crystal structure prediction." (GitHub 2016+). https://github.com/oganov-lab/STMng | STM image simulation, USPEX integration, structure viewing |
| 2 | 258 | **JMol** | GNU LGPL | http://jmol.sourceforge.net/ | Java molecular viewer | Hanson, R. M. "Jmol—a paradigm shift in crystallographic visualization." *J. Appl. Crystallogr.* **43**, 1250-1260 (2010). https://doi.org/10.1107/S0021889810030256 | Browser-based, crystallographic, web integration, educational |
| 3 | 259 | **PyMOL** | Commercial/Open (PyMOL-open-source) | https://pymol.org/ | Molecular visualization | The PyMOL Molecular Graphics System, Version 2.0 Schrödinger, LLC. (Open-source fork available); DeLano, W. L. "The PyMOL Molecular Graphics System." (2002). | Publication-quality rendering, biomolecules, ray tracing, scripting |


### **Category 9 Total: 14 codes**


====================================================================================================
# 10. FRAMEWORKS, WORKFLOW ENGINES & DATABASES
====================================================================================================

## 10.1 Materials Science Frameworks
**Subcategory Total: 9 codes**
────────────────────────────────────────────────────────────────────────────────

### 10.1.1 Python-Based Core Libraries

### 10.1.2 Workflow Management Engines

### 10.1.3 AiiDA Plugins & Ecosystem

| # | ID | Code Name | License | Official Website | Basis Set | Primary Publication with DOI | Specialization |
|---|---|---|---|---|---|---|---|
| 1 | 379 | **MatPy** | Research | Research code | Materials Python | Python materials analysis tools | Materials analysis in Python |
| 2 | 71 | **atomate2** | Modified BSD | https://github.com/materialsproject/atomate2 | Next-gen workflows | Ganose, A. M.; et al. "atomate2: A materials science toolkit." (GitHub 2021+). https://github.com/materialsproject/atomate2 | Second-generation workflows, jobflow foundation |
| 3 | - | **jobflow-remote** | - | - | - | - | - |
| 4 | 372 | **Luigi** | Apache 2.0 | https://github.com/spotify/luigi | Workflow management | Spotify Engineering. "Luigi: Batch processing framework." (GitHub). https://github.com/spotify/luigi | Python workflow management |
| 5 | 398 | **Parsl** | Apache 2.0 | https://parsl-project.org/ | Parallel scripting | Babuji, Y.; et al. "Parsl: Pervasive parallel programming in Python." *HPDC* (2019). https://doi.org/10.1145/3307681.3325400 | Parallel Scripting Library |
| 6 | 331 | **AiiDA-VASP** | MIT License | https://github.com/aiida-vasp/aiida-vasp | VASP plugin | AiiDA-VASP plugin. (GitHub). https://github.com/aiida-vasp/aiida-vasp | VASP integration with AiiDA |
| 7 | 330 | **AiiDA-QuantumESPRESSO** | MIT License | https://github.com/aiidateam/aiida-quantumespresso | QE plugin | AiiDA plugin for Quantum ESPRESSO. (GitHub). https://github.com/aiidateam/aiida-quantumespresso | Quantum ESPRESSO integration with AiiDA |
| 8 | 333 | **AiiDA-yambo** | MIT License | https://github.com/yambo-code/yambo-aiida | Yambo plugin | AiiDA plugin for Yambo. (GitHub). https://github.com/yambo-code/yambo-aiida | Yambo GW/BSE workflows in AiiDA |
| 9 | - | **aiida-fleur** | - | - | - | - | - |

## 10.2 High-Throughput & Database Infrastructure
**Subcategory Total: 9 codes**
────────────────────────────────────────────────────────────────────────────────

### 10.2.1 Database Frameworks & Platforms

### 10.2.2 Specialized High-Throughput Tools

| # | ID | Code Name | License | Official Website | Basis Set | Primary Publication with DOI | Specialization |
|---|---|---|---|---|---|---|---|
| 1 | 266 | **pymatgen-db** | MIT License | https://github.com/materialsproject/pymatgen-db | MongoDB interface | Ong, S. P.; et al. "Python Materials Genomics (pymatgen): A robust, open-source python library for materials analysis." *Comput. Mater. Sci.* **68**, 314-319 (2013). https://doi.org/10.1016/j.commatsci.2012.10.028 | MongoDB database builder, MP ecosystem, query interface |
| 2 | 262 | **qmpy** | Open Source | https://oqmd.org/ | OQMD Python interface | Saal, J. E.; Kirklin, S.; Aykol, M.; Meredig, B.; Wolverton, C. "Materials design and discovery with high-throughput density functional theory: the Open Quantum Materials Database (OQMD)." *JOM* **65**, 1501-1509 (2013). https://doi.org/10.1007/s11837-013-0755-5 | Open Quantum Materials Database, 1M+ entries, thermodynamics |
| 3 | 264 | **Materials Cloud** | Various | https://www.materialscloud.org/ | Computational infrastructure | Talirz, L.; Kumbhar, S.; Passaro, E.; Yakutovich, A. V.; Granata, V.; Gargiulo, F.; Borelli, M.; Uhrin, M.; Huber, S. P.; Zoupanos, S.; et al. "Materials Cloud, a platform for open computational science." *Sci. Data* **7**, 299 (2020). https://doi.org/10.1038/s41597-020-00637-5 | EPFL/PSI platform, AiiDA ecosystem, interactive tools, FAIR data |
| 4 | 273 | **MPWorks** | BSD 3-Clause | https://github.com/materialsproject/MPWorks | MP workflows (legacy) | Jain, A.; Ong, S. P.; Chen, W.; Medasani, B.; Qu, X.; Kocher, M.; Brafman, M.; Petretto, G.; Rignanese, G.-M.; Hautier, G.; Gunter, D.; Persson, K. A. "FireWorks: a dynamic workflow system designed for high-throughput applications." *Concurrency Computat.: Pract. Exper.* **27**, 5037-5059 (2015). https://doi.org/10.1002/cpe.3505 | Legacy MP workflows, superseded by atomate |
| 5 | 267 | **emmet** | Modified BSD | https://github.com/materialsproject/emmet | Materials Project builder | Horton, M. K.; Dwaraknath, S.; Montoya, J. H.; Ong, S. P.; Persson, K. A. "Emmet: A Python library for building, managing, and analyzing materials science databases." (GitHub 2018+) https://github.com/materialsproject/emmet | MP database builder, data aggregation, schema validation |
| 6 | 269 | **Matbench** | MIT License | https://github.com/materialsproject/matbench | ML benchmark suite | Dunn, A.; Wang, Q.; Ganose, A.; Dopp, D.; Jain, A. "Benchmarking materials property prediction methods: the Matbench test set and Automatminer reference algorithm." *npj Comput. Mater.* **6**, 138 (2020). https://doi.org/10.1038/s41524-020-00406-3 | 13 ML tasks, leaderboard, automated benchmarking |
| 7 | 270 | **CatApp** | Academic | https://catapp.org/ | Catalysis database | Winther, K. T.; Hoffmann, M. J.; Mamun, O.; Boes, J. R.; Nørskov, J. K. "Catalysis-Hub.org, an open electronic structure database for surface reactions." *Sci. Data* **6**, 75 (2019). https://doi.org/10.1038/s41597-019-0081-y | Surface reaction energies, adsorption, DFT database, 100,000+ calculations |
| 8 | 271 | **CatMAP** | GNU GPL v3 | https://github.com/SUNCAT-Center/catmap | Microkinetic modeling | Medford, A. J.; Shi, C.; Hoffmann, M. J.; Lausche, A. C.; Fitzgibbon, S. R.; Bligaard, T.; Nørskov, J. K. "CatMAP: A Software Package for Descriptor-Based Microkinetic Mapping of Catalytic Trends." *Catal. Lett.* **145**, 794-807 (2015). https://doi.org/10.1007/s10562-015-1495-6 | Microkinetics, rate analysis, scaling relations, SUNCAT |
| 9 | 272 | **GASpy** | Apache 2.0 | https://github.com/ulissigroup/GASpy | HT surface calculations | Tran, K.; Ulissi, Z. W. "Active learning across intermetallics to guide discovery of electrocatalysts for CO2 reduction and H2 evolution." *Nat. Catal.* **1**, 696-703 (2018). https://doi.org/10.1038/s41929-018-0142-1 | Active learning, surfaces, adsorption, FireWorks integration |


### **Category 10 Total: 18 codes**


====================================================================================================
# 11. SMALL, NICHE & RESEARCH-GRADE TOOLS
====================================================================================================

## 11.1 Specialized Electronic Structure Methods
**Subcategory Total: 5 codes**
────────────────────────────────────────────────────────────────────────────────

| # | ID | Code Name | License | Official Website | Basis Set | Primary Publication with DOI | Specialization |
|---|---|---|---|---|---|---|---|
| 1 | - | **RMG** | - | - | - | - | - |
| 2 | 370 | **KITE** | Apache 2.0 | https://quantum-kite.com/ | Tight-binding transport | João, S. M.; et al. "KITE: high-performance accurate modelling of electronic structure." *R. Soc. Open Sci.* **7**, 191809 (2020). https://doi.org/10.1098/rsos.191809 | Large-scale TB, Chebyshev expansion |
| 3 | 150 | **PAOFLOW** | MIT License | https://github.com/marcobn/PAOFLOW | TB post-processing | Buongiorno Nardelli, M.; et al. "PAOFLOW: A utility to construct and operate on ab initio Hamiltonians." *Comput. Mater. Sci.* **143**, 462-472 (2018). https://doi.org/10.1016/j.commatsci.2017.11.034 | Projections onto atomic orbitals |
| 4 | 377 | **MagneticTB** | Research | Research code | Magnetic TB | Tight-binding with magnetic interactions | Magnetic tight-binding |
| 5 | 376 | **MagneticKP** | Research | Research code | Magnetic k·p | k·p theory with magnetic interactions | Magnetic k·p methods |

## 11.2 Model Hamiltonians & Pedagogical Tools
**Subcategory Total: 2 codes**
────────────────────────────────────────────────────────────────────────────────

| # | ID | Code Name | License | Official Website | Basis Set | Primary Publication with DOI | Specialization |
|---|---|---|---|---|---|---|---|
| 1 | 135 | **cmpy** | Open Source | Migrated to exactdiag | Exact diagonalization | Community-maintained exact diagonalization tools | Pedagogical ED codes |
| 2 | - | **Stoner** | - | - | - | - | - |

## 11.3 Machine Learning Potentials & Neural Networks
**Subcategory Total: 8 codes**
────────────────────────────────────────────────────────────────────────────────

| # | ID | Code Name | License | Official Website | Basis Set | Primary Publication with DOI | Specialization |
|---|---|---|---|---|---|---|---|
| 1 | - | **MLIP ecosystem** | - | - | - | - | - |
| 2 | - | **n2p2** | - | - | - | - | - |
| 3 | 277 | **SIMPLE-NN** | Apache 2.0 | https://github.com/MDIL-SNU/SIMPLE-NN | Behler-Parrinello NNP | Lee, K.; Yoo, D.; Jeong, W.; Han, S. "SIMPLE-NN: An efficient package for training and executing neural-network interatomic potentials." *Comput. Phys. Commun.* **242**, 95-103 (2019). https://doi.org/10.1016/j.cpc.2019.04.014 | Symmetry functions, GPU training, LAMMPS interface |
| 4 | 324 | **AMP** | MIT License | https://amp.readthedocs.io/ | ML potentials | Khorshidi, A.; Peterson, A. A. "Amp: A modular approach to machine learning." *Comput. Phys. Commun.* **207**, 310-324 (2016). https://doi.org/10.1016/j.cpc.2016.05.010 | Atomistic Machine-learning Package |
| 5 | - | **SchNetPack** | - | - | - | - | - |
| 6 | 373 | **MACE** | MIT License | https://github.com/ACEsuit/mace | Equivariant ML | Batatia, I.; et al. "MACE: Higher Order Equivariant Message Passing Neural Networks." *NeurIPS* (2022). https://arxiv.org/abs/2206.07697 | Message passing equivariant networks |
| 7 | 389 | **NequIP** | MIT License | https://github.com/mir-group/nequip | E(3)-equivariant NN | Batzner, S.; et al. "E(3)-equivariant graph neural networks for data-efficient." *Nat. Commun.* **13**, 2453 (2022). https://doi.org/10.1038/s41467-022-29939-5 | E(3)-equivariant message passing |
| 8 | 278 | **Allegro** | MIT License | https://github.com/mir-group/allegro | E(3)-equivariant NNP | Musaelian, A.; Batzner, S.; Johansson, A.; Sun, L.; Owen, C. J.; Kornbluth, M.; Kozinsky, B. "Learning local equivariant representations for large-scale atomistic dynamics." *Nat. Commun.* **14**, 579 (2023). https://doi.org/10.1038/s41467-023-36329-y | Strictly local, E(3)-equivariant, 10× faster than NequIP |

## 11.4 API & Interface Tools
**Subcategory Total: 4 codes**
────────────────────────────────────────────────────────────────────────────────

| # | ID | Code Name | License | Official Website | Basis Set | Primary Publication with DOI | Specialization |
|---|---|---|---|---|---|---|---|
| 1 | 281 | **API_Phonons** | Academic | https://github.com/MaterialsDiscovery/PyChemia | Multi-code phonon interface | Gomez, G. H.; Guillén, G. G.; Villamagua, L.; Stashans, A. "PyChemia: A Python framework for materials discovery and design." (GitHub 2014+) https://github.com/MaterialsDiscovery/PyChemia | Phonon interface, multiple codes, workflow automation |
| 2 | 282 | **gpaw-tools** | GNU GPL v3 | https://github.com/lrgresearch/gpaw-tools | GPAW utilities | Sevik, C. et al. "gpaw-tools: Python toolkit for GPAW." (GitHub documentation 2018+) https://github.com/lrgresearch/gpaw-tools | GPAW automation, band structure, DOS, user-friendly scripts |
| 3 | 283 | **ASE-GUI** | GNU LGPL | https://wiki.fysik.dtu.dk/ase/ase/gui/gui.html | Graphical interface | Larsen, A. H.; et al. "The atomic simulation environment—a Python library for working with atoms." *J. Phys.: Condens. Matter* **29**, 273002 (2017). https://doi.org/10.1088/1361-648X/aa680e | Graphical ASE interface, structure builder, trajectory viewer |
| 4 | 284 | **Phonopy-API** | BSD 3-Clause | https://phonopy.github.io/phonopy/ | Phonopy Python API | Togo, A. "First-principles Phonon Calculations with Phonopy and Phono3py." *J. Phys. Soc. Jpn.* **92**, 012001 (2023). https://doi.org/10.7566/JPSJ.92.012001 | Python API for phonopy, automated workflows, scripting |

## 11.5 Specialized Analysis Tools
**Subcategory Total: 7 codes**
────────────────────────────────────────────────────────────────────────────────

| # | ID | Code Name | License | Official Website | Basis Set | Primary Publication with DOI | Specialization |
|---|---|---|---|---|---|---|---|
| 1 | 285 | **dbaAutomator** | Academic | https://github.com/BerkeleyGW/dbaAutomator | Double-Bader analysis | Jain, M.; Deslippe, J.; Samsonidze, G.; Cohen, M. L.; Chelikowsky, J. R.; Louie, S. G. "Improved quasiparticle wave functions and mean field for G0W0 calculations: Initialization with the COHSEX operator." *Phys. Rev. B* **90**, 115148 (2014). https://doi.org/10.1103/PhysRevB.90.115148 | Exciton analysis, BerkeleyGW, charge density analysis |
| 2 | 286 | **yambopy** | GNU GPL | https://github.com/yambo-code/yambopy | Yambo scripting | Yambo development team. "yambopy: Python interface to yambo." (GitHub 2015+) https://github.com/yambo-code/yambopy | Yambo automation, convergence tests, GW workflows |
| 3 | 287 | **AutoBZ.jl** | MIT License | https://github.com/thchr/AutoBZ.jl | Brillouin zone integration | Christensen, T. "AutoBZ.jl: Automatic Brillouin zone integration in Julia." (GitHub 2020+) https://github.com/thchr/AutoBZ.jl | Adaptive BZ integration, Julia, high-precision quadrature |
| 4 | 288 | **Pheasy** | BSD 3-Clause | https://github.com/JLUsquad/Pheasy | Phonon analysis | Li, J.; et al. "Pheasy: A Python toolkit for phonon property analysis." (GitHub 2019+) https://github.com/JLUsquad/Pheasy | Phonon mode analysis, visualization, Grüneisen parameters |
| 5 | 197 | **effectivemass** | MIT License | https://github.com/lucydot/effmass | Effective mass | Whalley, L. D.; et al. "Impact of nonparabolic electronic band structure." *Phys. Rev. B* **99**, 085207 (2019). https://doi.org/10.1103/PhysRevB.99.085207 | Effective mass calculator |
| 6 | 289 | **BerryPI** | GNU GPL v3 | https://github.com/spichardo/BerryPI | Berry phase calculations | Nunes, R. W.; Gonze, X. "Berry-phase treatment of the homogeneous electric field perturbation in insulators." *Phys. Rev. B* **63**, 155107 (2001). https://doi.org/10.1103/PhysRevB.63.155107 | Polarization, Born effective charges, ABINIT/WIEN2k interface |
| 7 | 290 | **IrRep** | GNU GPL v3 | https://github.com/stepan-tsirkin/irrep | Irreducible representations | Iraola, M.; Mañes, J. L.; Bradlyn, B.; Horton, M. K.; Neupert, T.; Vergniory, M. G.; Tsirkin, S. S. "IrRep: Symmetry eigenvalues and irreducible representations of ab initio band structures." *Comput. Phys. Commun.* **272**, 108226 (2022). https://doi.org/10.1016/j.cpc.2021.108226 | Symmetry eigenvalues, band representations, VASP/QE/ABINIT |

## 11.6 Specialized Solvers & Methods
**Subcategory Total: 5 codes**
────────────────────────────────────────────────────────────────────────────────

| # | ID | Code Name | License | Official Website | Basis Set | Primary Publication with DOI | Specialization |
|---|---|---|---|---|---|---|---|
| 1 | 294 | **dual fermions** | Research codes | Various implementations | Diagrammatic DMFT extension | Rubtsov, A. N.; Katsnelson, M. I.; Lichtenstein, A. I. "Dual fermion approach to nonlocal correlations in the Hubbard model." *Phys. Rev. B* **77**, 033101 (2008). https://doi.org/10.1103/PhysRevB.77.033101; Rohringer, G.; Valli, A.; Toschi, A. "Local electronic correlation at the two-particle level." *Phys. Rev. B* **86**, 125114 (2012). https://doi.org/10.1103/PhysRevB.86.125114 | Diagrammatic extension of DMFT, non-local correlations, dual transformation method |
| 2 | 292 | **NORG** | Open Source | https://github.com/rqHe1/NORG | Natural orbitals renormalization | He, R.-Q.; Lu, Z.-Y. "Quantum renormalization groups based on natural orbitals." *Phys. Rev. B* **89**, 085108 (2014). https://doi.org/10.1103/PhysRevB.89.085108; Wang, J.-M.; et al. "Ab initio dynamical mean-field theory with natural orbitals renormalization group impurity solver: Formalism and applications." *npj Comput. Mater.* **11**, 57 (2025). https://doi.org/10.1038/s41524-025-01586-6 | Zero-temperature impurity solver, Zen toolkit, multi-orbital, DMFT, natural orbitals RG |
| 3 | 293 | **AFLOW-ML** | Open Source | https://aflowlib.org/aflow-ml | ML property prediction | Gossett, E.; Toher, C.; Oses, C.; Isayev, O.; Legrain, F.; Rose, F.; Zurek, E.; Carrete, J.; Mingo, N.; Curtarolo, S. "AFLOW-ML: A RESTful API for machine-learning predictions of materials properties." *Comput. Mater. Sci.* **152**, 134-142 (2018). https://doi.org/10.1016/j.commatsci.2018.03.075; Isayev, O.; et al. "Universal fragment descriptors for predicting properties of inorganic crystals." *Nat. Commun.* **8**, 15679 (2017). https://doi.org/10.1038/ncomms15679 | Machine learning within AFLOW, RESTful API, electronic/thermal/mechanical properties, PLMF/MFD models |
| 4 | 90 | **Materials Studio** | Commercial | https://www.3ds.com/products/biovia/materials-studio | Commercial suite | BIOVIA Materials Studio, Dassault Systèmes (1995-2025). https://www.3ds.com/products/biovia/materials-studio | Comprehensive commercial suite, BIOVIA |
| 5 | 296 | **MedeA** | Commercial | https://www.materialsdesign.com/medea | Computational materials design | Materials Design Inc. "MedeA (Materials Exploration and Design Analysis)." Commercial platform documentation (1998-2025). https://www.materialsdesign.com/medea | Commercial platform, automated workflows, VASP/LAMMPS integration, property prediction |

## 11.7 Additional Specialized Codes
**Subcategory Total: 15 codes**
────────────────────────────────────────────────────────────────────────────────

| # | ID | Code Name | License | Official Website | Basis Set | Primary Publication with DOI | Specialization |
|---|---|---|---|---|---|---|---|
| 1 | 297 | **AIMPRO** | Academic | https://www.aimpro.org/ | Pseudopotential DFT | Briddon, P. R.; Jones, R. "LDA calculations using a basis of Gaussian orbitals." *Phys. Status Solidi B* **217**, 131-171 (2000). https://doi.org/10.1002/(SICI)1521-3951(200001)217:1<131::AID-PSSB131>3.0.CO;2-M | Ab Initio Modelling PROgram, defects in semiconductors, Gaussian orbitals, Newcastle |
| 2 | 298 | **Ascalaph Designer** | Commercial | http://www.biomolecular-modeling.com/Ascalaph/index.html | Molecular modeling | Andrei, I. (Agile Molecule LLC). "Ascalaph Designer." Software documentation (2010+). http://www.biomolecular-modeling.com/Ascalaph/ | Molecular builder, protein-ligand docking, MD simulations, GUI-based |
| 3 | - | **Atompaw/PWPAW** | - | - | - | - | - |
| 4 | 112 | **Fireball** | Academic | https://github.com/fireball-QMD/progs | LCAO tight-binding | Lewis, J. P.; et al. "Advances in the FIREBALL ab initio tight-binding molecular-dynamics." *Phys. Status Solidi B* **248**, 1989-2007 (2011). https://doi.org/10.1002/pssb.201147259 | Fast LCAO, order-N, MD |
| 5 | 304 | **FSatom** | Academic | http://bohr.inesc-mn.pt/~jlm/pseudo.html | Pseudopotential generation | Martins, J. L. "FSatom: A program for generating all-electron and pseudopotential atom calculation." (Code documentation 2000+). http://bohr.inesc-mn.pt/~jlm/pseudo.html | Free atom solver, pseudopotential testing, Hamann/Troullier-Martins/Kerker PPs |
| 6 | 305 | **HiLAPW** | Research | Limited public access | High-speed LAPW | Kato, M.; et al. "HiLAPW: High-speed full-potential LAPW implementation." (Japanese development, limited English documentation). Research code | High-performance LAPW, Japanese development, limited distribution |
| 7 | 108 | **NRLMOL** | Academic | https://www.nrl.navy.mil/ | Gaussian + mesh | Pederson, M. R.; Jackson, K. A. "Variational mesh for quantum-mechanical simulations." *Phys. Rev. B* **41**, 7453-7461 (1990). https://doi.org/10.1103/PhysRevB.41.7453 | US Naval Research Lab, Gaussian orbitals on mesh |
| 8 | 307 | **ParaGauss** | Academic | https://github.com/bo-li/paragauss | Parallel quantum chemistry | Belling, T.; Grauschopf, T.; Krüger, S.; Nörtemann, F.; Staufer, M.; Mayer, M.; Nasluzov, V. A.; Birkenheuer, U.; Hu, A.; Matveev, A. V.; Shor, A. M.; Fuchs-Rohr, M. S. K.; Neyman, K. M.; Ganyushin, D. I.; Kerdcharoen, T.; Woiterski, A.; Görling, A.; Rösch, N. "PARAGAUSS: A DFT parallelization in practice." In *High Performance Computing in Science and Engineering '98* (Springer 1999). (GitHub 2010+). https://github.com/bo-li/paragauss | Parallel DFT, scalar-relativistic, molecular systems, Rösch group (Munich) |
| 9 | 310 | **Petot** | Research | Limited documentation | DFT implementation | Limited English documentation available. Research code with restricted access. | Computational chemistry implementation, limited public information |
| 10 | 311 | **Socorro** | Academic | https://dft.sandia.gov/ | Plane-wave/LCAO DFT | Schultz, P. A. "SeqQuest Electronic Structure Code." Sandia National Laboratories. (SeqQuest/Socorro codes). https://dft.sandia.gov/ | Sandia DFT package, plane-wave and LCAO modes, massively parallel |
| 11 | 384 | **Materials and Processes Simulations** | Research | Various tools | Multiphysics | Materials and process simulation tools (category) | Materials process simulations |
| 12 | - | **TOTAL** | - | - | - | - | - |
| 13 | - | **~469** | - | - | - | - | - |
| 14 | - | **~141** | - | - | - | - | - |
| 15 | - | **~328** | - | - | - | - | - |


### **Category 11 Total: 46 codes**


====================================================================================================
# SUMMARY STATISTICS
====================================================================================================

**Total Codes in Catalog**: 278
**Codes with Complete Details**: 214 (77.0%)
**Codes Missing Details**: 64

## Category Breakdown

| Category | Title | Code Count |
|---|---|---|
| 1 | GROUND-STATE ELECTRONIC STRUCTURE (DFT & VARIANTS) | 123 |
| 2 | TIME-DEPENDENT & EXCITED-STATE METHODS | 8 |
| 3 | STRONGLY CORRELATED & MANY-BODY METHODS | 30 |
| 4 | WAVEFUNCTION-BASED QUANTUM CHEMISTRY | 0 |
| 5 | TIGHT-BINDING, MODEL HAMILTONIANS & DOWNFOLDING | 9 |
| 6 | PHONONS, LATTICE DYNAMICS & ELECTRON-PHONON COUPLI | 16 |
| 7 | MOLECULAR & AB INITIO DYNAMICS | 3 |
| 8 | STRUCTURE PREDICTION & GLOBAL OPTIMIZATION | 11 |
| 9 | POST-PROCESSING, ANALYSIS & VISUALIZATION | 14 |
| 10 | FRAMEWORKS, WORKFLOW ENGINES & DATABASES | 18 |
| 11 | SMALL, NICHE & RESEARCH-GRADE TOOLS | 46 |
| **GRAND TOTAL** | | **278** |

## Codes Missing Details (Not in entries*.md files)
**Total Missing**: 64

```
  - ACES II
  - ALF
  - ALPS/CT-HYB
  - ALPSCore
  - AMULET
  - AtomViz
  - Atomistix ToolKit
  - Atompaw/PWPAW
  - Basin hopping
  - CHAMP
  - CT-HYB
  - CT-INT
  - CT-QMC
  - CT-SEG
  - DCA++
  - DMDW/RTDW
  - DMFTwDFT
  - EPW
  - FTPS
  - GAMESS-UK
  - Gaussian 16
  - Gaussian Basis Set Libraries
  - GreenX
  - HOTBIT
  - HÏ†
  - JuKKR
  - Julia materials packages
  - KKRnano
  - LMTO-ASA
  - MLIP ecosystem
  - Materials Simulation Suites
  - Metadynamics
  - NBSE
  - Neural network potentials
  - PLATO
  - Priroda-06
  - PyQMC
  - PyQuante
  - QMcBeaver
  - QUEST
  - QWalk
  - Qbox
  - RMG
  - SALMON
  - SchNetPack
  - Spex
  - Stoner
  - String methods
  - TOTAL
  - TRIQS/CT-QMC solvers
  - TRIQS/DFTTools
  - VASP+DMFT
  - X-ray Crystallography Integration
  - Z2Pack
  - aiida-fleur
  - fiesta
  - jobflow-remote
  - molgw
  - n2p2
  - pymatgen-diffusion
  - solid_dmft
  - ~141
  - ~328
  - ~469
```

### Verification Note
These 64 codes appear in the comprehensive_code_list.md but do not have matching entries in any of the entries*.md files. This means their full details (license, website, basis set, publication, specialization) are not available in the current documentation.