## Loc_Aper.f90 – Localization and Input Generation Module for Local Post-HF Treatments

This module was developed as part of my master’s thesis at the Theoretical Chemistry Department of Humboldt University of Berlin, under the supervision of PD Dr. Denis Usvyat. It is implemented within the Cryscor program and is designed to enable local post-Hartree–Fock treatments of an aperiodic fragment.

### Features

The module sets up the local approximations required for the correlation treatment of a predefined fragment, including:

- `my_loc`: Constructs the localized molecular orbitals (LMOs) for the occupied space using the Pipek–Mezey procedure.
- `Frag_PAOs`: Constructs the projected atomic orbitals (PAOs) used to represent the virtual space and normalizes them in the orthonormal PAO basis.
- `Frag_domains`:  Builds minimal orbital domains for each localized orbital and extends them based on bond connectivity (number of bonds) and/or interatomic distances. These domains define the local environment used for correlation.
- `Frag_pair_list`: Constructs the pair list by classifying orbital pairs as strong, weak, or distant, based on the minimum distance between atoms belonging to their respective orbital domains.
- `Frag_Fock`: Constructs the transformation matrices that project the Fock matrix from the atomic orbital (AO) basis into the local basis. The occupied part is transformed to the LMO basis and the virtual part to the PAO basis. The routine also removes linear dependencies in the PAO space.
- `write_frag_molpro_input`: Prepares the input file required to run local post-HF calculations within the Local Integrated Tensor Framework (LITF), part of the Molpro package. The file includes all localized quantities for the fragment: orbital domains, pair list, Fock matrices, overlap matrices, transformation matrices, and electron repulsion integrals.

---

Author: Lily Barta  
Email: lily.barta@etu.u-bordeaux.fr  
International Master’s Program in Physical Chemistry and Chemical Physics (University of Bordeaux)
