# Adaptive-SCF-ELK-7.2.42
An adaptive preconditioning scheme that can identify the long-wavelength divergence behavior of the Jacobian and choose the appropriate preconditioning method during the SCF iteration.\\
The adaptive preconditioning algorithm is implemented based on the Elk-7.2.42 package.

This algorithm addes the following programs (in package):\\
msbessel1.f90 	msbessel2.f90	kgenrik.f90	kpotcoul.f90 
kzpotcoul.f90	kgenzvclmt.f90	kzpotclmt.f90	mix_final.f90  
mix_init.f90  	mix_scf.f90	pmix_scf.f90	modplmix.f90

In addition, we also modify these files of Elk-7.2.42 package:  readinput.f90 and gndstate.f90.

msbessel1.f90, msbessel2.f90, kgenrik.f90, kpotcoul.f90, kzpotcoul.f90, kgenzvclmt.f90, and kzpotclmt.f90 are related to the codes mentioned in the manuscript for solving the screened Poisson equation; mix_final.f90, mix_init.f90, mix_scf.f90, pmix_scf.f90, and modplmix.f90 are related to the codes mentioned in the manuscript for updating the subspace and calculation of posterior indicator based on Pulay mixing. 

The above files need to be placed in the elk-7.2.42/src directory.

Regarding the compilation of the program, we provide the modified elk-7.2.42/make.inc and elk-7.2.42/src/Makefile as references.
 
We use the LAPACK and BLAS libraries of Intel in make.inc (MKLROOT, MKL_PATH, and LIB_LPK). It is worth noting that we employed Intel compilers, and users need to modify the corresponding path according to their compiler.

For the adaptive preconditioning algorithm, we added two parameters: klambda and premix.

klambda (real) 
Default value: klambda = 0.6

In the screened Poisson equation $(\nabla^2 - \lambda^2)V(\textbf{r})= - 4\pi\rho(\textbf{r})$, $ klambda = \lambda $. 
klambda is usually set to Thomas-Fermi wave number. In principle，it represents the inverse of a typical length scale over which an individual charged particle exerts a notable effect.  klambda = 0.6 is a common setting.

premix (integer)
Default value: premix = 1

premix is the type of preconditioning:
0	:  non-preconditioning;
1	:  adaptive preconditioning;
2	:  Kerker preconditioning.
