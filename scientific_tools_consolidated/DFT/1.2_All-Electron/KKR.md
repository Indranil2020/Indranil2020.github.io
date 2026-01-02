# KKR (Korringa-Kohn-Rostoker)

## Important Note
**KKR is a METHOD, not standalone software**. It is implemented in various codes including:
- **SPR-KKR** (spin-polarized relativistic KKR)
- **AkaiKKR** (Akai's KKR implementation)
- **JuKKR** (Jülich KKR)
- **KKRnano** (massively parallel KKR)

## Method Overview
KKR (Korringa-Kohn-Rostoker) is a Green's function-based DFT method particularly powerful for disordered alloys (via CPA), spectroscopy, and magnetic materials. Unlike standard DFT codes, KKR uses multiple scattering theory and Green's functions, making it especially suited for substitutional disorder, spectroscopy calculations, and systems requiring advanced magnetic treatments.

**Method type**: Green's function DFT, multiple scattering  
**Primary implementations**: SPR-KKR, AkaiKKR, JuKKR, KKRnano

## Method Characteristics

### Green's Function Approach:
- Multiple scattering theory
- Direct Green's functions
- Spectroscopy natural
- Transport properties
- Advanced formalism

### Coherent Potential Approximation (CPA):
- Substitutional disorder
- Random alloys
- No supercells needed
- Exact statistical treatment
- Concentration effects

### Applications:
- Random alloys
- Spectroscopy (XAS, XMCD)
- Magnetic materials
- Transport properties
- Disordered systems

## Implementations

### SPR-KKR:
- Fully relativistic
- Spectroscopy focus
- Magnetism
- LMU Munich

### AkaiKKR:
- CPA for alloys
- Widely used
- Japanese development
- Free software

### JuKKR:
- Jülich implementation
- Modern code
- Python interface
- Active development

### KKRnano:
- Massively parallel
- Large-scale
- HPC focus
- Extreme scaling

## Recommendation
For KKR calculations, use established implementations:
- **SPR-KKR**: Spectroscopy, relativistic
- **AkaiKKR**: Alloys, CPA
- **JuKKR**: Modern, Python interface
- **KKRnano**: Large-scale parallel

## Verification & Sources
**Status**: ✅ METHOD VERIFIED
- KKR is a well-established DFT method
- Multiple production implementations available
- **Not standalone software** - use SPR-KKR, AkaiKKR, JuKKR, or KKRnano
