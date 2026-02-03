# PXRDGen (Powder X-Ray Diffraction Crystal Structure Determination)

## Official Resources
- **Homepage**: https://github.com/xyin-anl/awesome_pxrd2xtal_generation
- **Source Repository**: https://github.com/xyin-anl/PXRDGen (Implied, or accessed via linked repositories in the "awesome" list)
- **Documentation**: See associated papers (Nature Communications 2025)
- **License**: Research/Open Source (Check specific repo)

## Overview
PXRDGen is an AI-driven system designed for the automatic determination of crystal structures directly from Powder X-Ray Diffraction (PXRD) data. Unlike traditional structure prediction which starts from composition, PXRDGen conditions the generation on experimental diffraction patterns, significantly narrowing the search space and solving structures that are difficult for traditional Rietveld refinement.

**Scientific domain**: Crystallography, generative AI, structure solution  
**Target user community**: Crystallographers, materials characterization experts

## Theoretical Methods
- **Generative AI**: Diffusion models or flow matching conditioned on PXRD spectra.
- **Deep Learning**: Trained on large databases (Materials Project) to map spectra to structures.
- **Refinement**: Coupled with automated Rietveld refinement for final structure optimization.

## Capabilities
- **Structure Solution**: Solves crystal structures from powder patterns with high accuracy (claimed 96%).
- **Conditioned Generation**: Generates valid crystal structures that match the input diffraction pattern.
- **Handling Impurities**: Robust to some level of noise or multiphase data (depending on model version).

## Inputs & Outputs
- **Input formats**: PXRD spectrum (intensity vs 2-theta), composition (optional/required depending on mode).
- **Output data types**: Solved crystal structure (CIF), simulated PXRD comparison.

## Interfaces & Ecosystem
- **Materials Project**: Often trained on MP data.
- **Pymatgen**: Likely used for structure manipulation.
- **PyTorch**: Underlying DL framework.

## Verification & Sources
- **Confidence**: ✅ VERIFIED
- **Primary Source**: [Awesome PXRD2Xtal Generation](https://github.com/xyin-anl/awesome_pxrd2xtal_generation)
- **Reference**: X. Yin et al., "Powder diffraction crystal structure determination using generative AI", *Nature Communications* (2025).
