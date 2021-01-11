# CF_Hamiltonian
A crystal field module for GRASP developed by G. Gaigalas (Vilnius University, Lithuania) and D. Kato (National Institute for Fusion Science, Japan).

## Introduction
The latest version of the GRASP2018 package [Froese Fischer et al. (2019)], based on the multiconfigurational Dirac–Hartree–Fock method, is extended to account for effects of crystal fields in complex systems. Instead of using the simplified treatment of the crystal field effects based on the Stevens’ operator-equivalent method the program uses the fully ab-initio method in which the external ions are treated as point charges at fixed positions. In addition, examples of how to use the CF_Hamiltonian program are given in source directory grasp2018/src/appl/CF_Hamiltonian/Sample_Runs.

## Program summary
**Program Title:** CF_Hamiltonian  
**CPC Library link to program files:** https://doi.org/10.17632/fksxwwjbx6.1  
**Licensing provisions:** MIT license  
**Programming language:** Fortran 95  
**External routines/libraries used:** Grasp2018 modules: `Libmod`, `Lib9290`, `Librang90`; Grasp2018 routines: `starttime`, `setdbg`, `getmixblock`, `getmixa`, `getmixc`, `setmc`, `factt`, `setcon`, `setcsla`, `stoptime`; and the `Lapack` library  
**Nature of problem:** The CF_Hamiltonian program is designed as a part of the Grasp2018 package for the computation of Stark splitting in crystal field in the point charge crystal field approximation.  
**Solution method:** The point charge crystal field approach is used. It allows user to include different Atomic State Functions (ASF) mixing such as ASF mix with the same total J values, ASF mixing with different total J values, ASF mixing with different parities.  
**Additional comments including restrictions and unusual features:** The restrictions of the program are coming from the restrictions of Grasp2018 package and it is suitable for systems for which the point charge crystal field approximation is appropriate. The Stark level splitting of the atomic energy terms in the point-charge crystal field approach is performed by the program `CF_Hamiltonian`.  

Further details are given in the CPC article included under the `docs/` directory. Note that the installation procedure of the development version of GRASP may differ from how it is decribed in the paper, which is based on the CPC 2019 version of GRASP (GRASP2018).
