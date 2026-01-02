#!/usr/bin/env python3
"""
Generate skeleton .md files for all 372 tools
Zero-hallucination protocol: Create structure only, populate from verified sources
"""

import os
from pathlib import Path

# Tool list with category mapping
tools_data = {
    "DFT": [
        ("CASTEP", "CONFIRMED", "https://www.castep.org/"),
        ("CP2K", "CONFIRMED", "https://www.cp2k.org/"),
        ("CPMD", "CONFIRMED", "https://www.cpmd.org/"),
        ("GPAW", "CONFIRMED", "https://wiki.fysik.dtu.dk/gpaw/"),
        ("JDFTx", "CONFIRMED", "https://jdftx.org/"),
        ("Qbox", "CONFIRMED", "http://qboxcode.org/"),
        ("PARSEC", "VERIFIED", "http://parsec.ices.utexas.edu/"),
        ("PARATEC", "LOW_CONF", "UNKNOWN"),
        ("SPARC", "LOW_CONF", "https://sparc-x.github.io/"),
        ("RMGDFT", "LOW_CONF", "https://github.com/RMGDFT/rmgdft"),
        ("ABACUS", "LOW_CONF", "https://github.com/deepmodeling/abacus-develop"),
        ("PWPAW", "LOW_CONF", "UNKNOWN"),
        ("TBPW", "LOW_CONF", "UNKNOWN"),
        ("PROFESS", "LOW_CONF", "UNKNOWN"),
        ("MADNESS", "LOW_CONF", "UNKNOWN"),
        ("OpenAtom", "LOW_CONF", "UNKNOWN"),
        ("PWDFT", "LOW_CONF", "UNKNOWN"),
        ("PLATO", "LOW_CONF", "UNKNOWN"),
        ("NESSIE", "VERIFIED", "UNKNOWN"),
        ("DFT-FE", "LOW_CONF", "UNKNOWN"),
        ("WIEN2k", "CONFIRMED", "https://www.wien2k.at/"),
        ("Elk", "CONFIRMED", "http://elk.sourceforge.net/"),
        ("Fleur", "CONFIRMED", "https://www.flapw.de/"),
        ("exciting", "CONFIRMED", "https://exciting-code.org/"),
        ("Questaal", "CONFIRMED", "https://questaal.org/"),
        ("RSPt", "VERIFIED", "UNKNOWN"),
        ("SPR-KKR", "VERIFIED", "UNKNOWN"),
        ("JuKKR", "VERIFIED", "https://github.com/JuDFTteam/JuKKR"),
        ("KKRnano", "LOW_CONF", "UNKNOWN"),
        ("AkaiKKR", "LOW_CONF", "UNKNOWN"),
        ("LMTO-ASA", "LOW_CONF", "UNKNOWN"),
        ("FPLO", "LOW_CONF", "https://www.fplo.de/"),
        ("KKR", "LOW_CONF", "UNKNOWN"),
        ("LMTO", "LOW_CONF", "UNKNOWN"),
        ("FHI-aims", "CONFIRMED", "https://fhi-aims.org/"),
        ("SIESTA", "CONFIRMED", "https://siesta-project.org/siesta/"),
        ("OpenMX", "CONFIRMED", "http://www.openmx-square.org/"),
        ("CONQUEST", "CONFIRMED", "https://www.order-n.org/"),
        ("ONETEP", "CONFIRMED", "https://onetep.org/"),
        ("BigDFT", "CONFIRMED", "https://bigdft.org/"),
        ("CRYSTAL", "CONFIRMED", "http://www.crystal.unito.it/"),
        ("ADF", "VERIFIED", "https://www.scm.com/amsterdam-modeling-suite/adf/"),
        ("DMol3", "VERIFIED", "UNKNOWN"),
        ("deMon2k", "LOW_CONF", "UNKNOWN"),
        ("RMG", "LOW_CONF", "UNKNOWN"),
        ("ORCA", "CONFIRMED", "https://orcaforum.kofo.mpg.de/"),
        ("Gaussian", "CONFIRMED", "https://gaussian.com/"),
        ("PySCF", "CONFIRMED", "https://pyscf.org/"),
        ("PSI4", "CONFIRMED", "https://psicode.org/"),
        ("Molpro", "CONFIRMED", "https://www.molpro.net/"),
        ("NWChem", "CONFIRMED", "http://www.nwchem-sw.org/"),
        ("TURBOMOLE", "CONFIRMED", "UNKNOWN"),
        ("Q-Chem", "CONFIRMED", "https://www.q-chem.com/"),
        ("GAMESS", "VERIFIED", "https://www.msg.chem.iastate.edu/gamess/"),
        ("Dalton", "VERIFIED", "https://www.daltonprogram.org/"),
        ("DIRAC", "VERIFIED", "https://diracprogram.org/"),
        ("CFOUR", "CONFIRMED", "https://www.cfour.de/"),
        ("MRCC", "CONFIRMED", "https://www.mrcc.hu/"),
        ("OpenMolcas", "CONFIRMED", "https://www.molcas.org/"),
        ("BAGEL", "CONFIRMED", "https://nubakery.org/"),
        ("Columbus", "VERIFIED", "UNKNOWN"),
        ("ACES", "LOW_CONF", "http://www.qtp.ufl.edu/ACES/"),
        ("ExaChem", "LOW_CONF", "https://github.com/exachem"),
        ("Quantum-Package", "LOW_CONF", "https://quantumpackage.github.io/qp2/"),
        ("CheMPS2", "LOW_CONF", "https://github.com/SebWouters/CheMPS2"),
        ("SlowQuant", "LOW_CONF", "https://github.com/slowquant/slowquant"),
        ("BDF", "LOW_CONF", "UNKNOWN"),
        ("eT", "LOW_CONF", "UNKNOWN"),
        ("CC4S", "UNCERTAIN", "UNKNOWN"),
        ("ACES-III", "UNCERTAIN", "UNKNOWN"),
        ("Molcas", "UNCERTAIN", "UNKNOWN"),
        ("DFTB+", "CONFIRMED", "https://www.dftbplus.org/"),
        ("xTB", "CONFIRMED", "https://github.com/grimme-lab/xtb"),
        ("HOTBIT", "LOW_CONF", "UNKNOWN"),
        ("MOPAC", "LOW_CONF", "UNKNOWN"),
        ("AMS-DFTB", "LOW_CONF", "https://www.scm.com/"),
        ("Materials-Studio", "LOW_CONF", "UNKNOWN"),
        ("Medea", "LOW_CONF", "UNKNOWN"),
        ("FLAPW", "LOW_CONF", "UNKNOWN"),
        ("FlapwMBPT", "LOW_CONF", "UNKNOWN"),
        ("DFT-F", "LOW_CONF", "UNKNOWN"),
    ],
    "TDDFT": [
        ("Octopus", "CONFIRMED", "https://octopus-code.org/"),
        ("SALMON", "CONFIRMED", "https://salmon-tddft.jp/"),
        ("Yambo", "CONFIRMED", "http://www.yambo-code.org/"),
        ("turboTDDFT", "VERIFIED", "UNKNOWN"),
        ("PyTDDFT", "LOW_CONF", "UNKNOWN"),
        ("TDAP", "UNCERTAIN", "UNKNOWN"),
        ("WEST", "CONFIRMED", "https://west-code.org/"),
        ("Spex", "CONFIRMED", "UNKNOWN"),
        ("SternheimerGW", "LOW_CONF", "UNKNOWN"),
        ("Fiesta", "VERIFIED", "UNKNOWN"),
        ("molgw", "LOW_CONF", "UNKNOWN"),
        ("GreenX", "LOW_CONF", "UNKNOWN"),
        ("SAX", "LOW_CONF", "UNKNOWN"),
        ("OCEAN", "VERIFIED", "UNKNOWN"),
        ("NBSE", "LOW_CONF", "UNKNOWN"),
        ("DP-Code", "LOW_CONF", "UNKNOWN"),
        ("DP-4", "LOW_CONF", "UNKNOWN"),
        ("pyGWBSE", "LOW_CONF", "UNKNOWN"),
    ],
    "DMFT": [
        ("TRIQS-DFTTools", "CONFIRMED", "https://triqs.github.io/dft_tools/"),
        ("TRIQS-cthyb", "VERIFIED", "UNKNOWN"),
        ("solid_dmft", "LOW_CONF", "UNKNOWN"),
        ("w2dynamics", "CONFIRMED", "https://github.com/w2dynamics/w2dynamics"),
        ("DCore", "CONFIRMED", "https://issp-center-dev.github.io/DCore/"),
        ("iQIST", "CONFIRMED", "UNKNOWN"),
        ("EDMFTF", "CONFIRMED", "https://hauleweb.rutgers.edu/"),
        ("ComDMFT", "CONFIRMED", "https://github.com/ComDMFT/ComDMFT"),
        ("ComCTQMC", "VERIFIED", "UNKNOWN"),
        ("ComRISB", "LOW_CONF", "UNKNOWN"),
        ("DMFTwDFT", "VERIFIED", "https://github.com/DMFTwDFT-project/DMFTwDFT"),
        ("AMULET", "LOW_CONF", "UNKNOWN"),
        ("Rutgers-DMFT", "LOW_CONF", "https://www.physics.rutgers.edu/~haule/CODES/"),
        ("ALPS", "VERIFIED", "UNKNOWN"),
        ("ALPSCore", "LOW_CONF", "UNKNOWN"),
        ("GTM", "UNCERTAIN", "UNKNOWN"),
        ("NRGLjubljana", "LOW_CONF", "UNKNOWN"),
        ("opendf", "LOW_CONF", "https://github.com/CQMP/opendf"),
        ("Kondo", "UNCERTAIN", "UNKNOWN"),
        ("COMSUITE", "LOW_CONF", "UNKNOWN"),
        ("CT-HYB", "VERIFIED", "UNKNOWN"),
        ("CT-QMC", "VERIFIED", "UNKNOWN"),
        ("CT-INT", "LOW_CONF", "UNKNOWN"),
        ("CT-SEG", "LOW_CONF", "UNKNOWN"),
        ("HubbardPhi", "LOW_CONF", "UNKNOWN"),
        ("EDIpack", "LOW_CONF", "UNKNOWN"),
        ("FTPS", "LOW_CONF", "UNKNOWN"),
        ("Pomerol", "VERIFIED", "UNKNOWN"),
        ("ITensor", "LOW_CONF", "UNKNOWN"),
        ("TeNPy", "LOW_CONF", "UNKNOWN"),
        ("Block", "LOW_CONF", "UNKNOWN"),
        ("DMRG++", "UNCERTAIN", "UNKNOWN"),
        ("NORG", "LOW_CONF", "UNKNOWN"),
        ("Dual-fermions", "LOW_CONF", "UNKNOWN"),
        ("EDRIXS", "LOW_CONF", "UNKNOWN"),
        ("exactdiag", "LOW_CONF", "UNKNOWN"),
        ("HubbardFermiMatsubara", "LOW_CONF", "UNKNOWN"),
    ],
    "QMC": [
        ("CASINO", "CONFIRMED", "https://vallico.net/casinoqmc/"),
        ("TurboRVB", "CONFIRMED", "https://github.com/sissaschool/turborvb"),
        ("ALF", "CONFIRMED", "https://alf.physik.uni-wuerzburg.de/"),
        ("CHAMP", "VERIFIED", "UNKNOWN"),
        ("QWalk", "VERIFIED", "http://www.qwalk.org/"),
        ("PyQMC", "LOW_CONF", "UNKNOWN"),
        ("QMcBeaver", "LOW_CONF", "UNKNOWN"),
        ("QUEST", "LOW_CONF", "https://github.com/QUEST-QMC/QUEST"),
        ("DCA++", "LOW_CONF", "UNKNOWN"),
        ("NECI", "LOW_CONF", "UNKNOWN"),
        ("HANDE", "LOW_CONF", "UNKNOWN"),
        ("ph-AFQMC", "LOW_CONF", "https://github.com/ph-AFQMC/ph-AFQMC"),
        ("qmclib", "UNCERTAIN", "UNKNOWN"),
        ("ZTC", "UNCERTAIN", "UNKNOWN"),
    ],
}

def create_skeleton_md(category, tool_name, confidence, resource):
    """Create skeleton markdown file for a tool"""
    
    template = f"""# {tool_name}

## Official Resources
- Homepage: {resource if resource != "UNKNOWN" else "UNKNOWN - Requires verification"}
- Documentation: UNKNOWN - Requires verification
- Source Repository: UNKNOWN - Requires verification
- License: UNKNOWN - Requires verification

## Overview
**Confidence Level**: {confidence}
**Status**: Documentation pending

[TO BE COMPLETED]

## Theoretical Methods
[TO BE COMPLETED - Requires verification from official sources]

## Capabilities (CRITICAL)
[TO BE COMPLETED - Only verified capabilities from official documentation]

**Sources**: Pending verification

## Inputs & Outputs
**Input formats**: [TO BE COMPLETED]

**Output data types**: [TO BE COMPLETED]

## Interfaces & Ecosystem
[TO BE COMPLETED - Requires verification]

## Limitations & Known Constraints
[TO BE COMPLETED - Requires official documentation review]

## Verification & Sources
**Primary sources**: [TO BE VERIFIED]

**Secondary sources**: [TO BE VERIFIED]

**Confidence**: {confidence}

**Verification status**: ⏸️ PENDING
- Official homepage: {'UNKNOWN' if resource == 'UNKNOWN' else 'TO BE VERIFIED'}
- Documentation: TO BE VERIFIED
- Capabilities: TO BE VERIFIED
"""
    
    return template

def main():
    base_dir = Path("/home/niel/git/Indranil2020.github.io/scientific_tools_consolidated")
    
    created_count = 0
    skipped_count = 0
    
    for category, tools in tools_data.items():
        category_dir = base_dir / category
        category_dir.mkdir(exist_ok=True)
        
        for tool_name, confidence, resource in tools:
            # Clean filename
            filename = tool_name.replace("/", "-").replace(" ", "-") + ".md"
            filepath = category_dir / filename
            
            # Skip if already exists
            if filepath.exists():
                skipped_count += 1
                print(f"SKIP: {filepath} (already exists)")
                continue
            
            # Create skeleton
            content = create_skeleton_md(category, tool_name, confidence, resource)
            filepath.write_text(content)
            created_count += 1
            print(f"CREATE: {filepath}")
    
    print(f"\n=== Summary ===")
    print(f"Created: {created_count} files")
    print(f"Skipped: {skipped_count} files (already exist)")
    print(f"Total: {created_count + skipped_count} files")

if __name__ == "__main__":
    main()
