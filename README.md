# rgmresir
MATLAB codes for performing gmres-based iterative refinement with recycling

This code can be used to reproduce the experiments in 

GCRO_DR codes are adapted from the codes developed in Parks, Michael L., et al. "Recycling Krylov subspaces for sequences of linear systems.", SIAM Journal on Scientific Computing 28.5 (2006): 1651-1674. The library and associated functions of GCRO-DR are available at https://www.sandia.gov/-mlparks/software/.
 
## Included MATLAB files
* **_chop.m, float_params.m, lutx_chop.m,  trisol.m, and roundit.m_** are functions in chop library that simulate half precision. The library and associated functions are available at https://github.com/higham/chop and https://github.com/SrikaraPranesh/LowPrecision\_Simulation.

* **_scale_diag_2side.m_** is a function that performs two-sided diagonal scaling of matrix. It is available at https://github.com/SrikaraPranesh/Multi_precision_NLA_kernels.

* **_lp_matvec.m_** is a function that performs matvec in low precision. It is available at https://github.com/SrikaraPranesh/Multi_precision_NLA_kernels.

* **_rgmresir.m_** is a function that performs GMRES-based iterative refinement with recycling in three precisions.

* **_gmresir.m_** is a function that performs GMRES-based iterative refinement in three precisions with an extra higher precision in factorization and preconditioner steps of GMRES.

* **_gmres_sd.m and gmres_dq.m_** are functions that run left-preconditioned GMRES using precisions single/double and double/quad. Application of the preconditioned coefficient matrix to a vector and the preconditioner to the right-hand-side vector are performed in the higher precision; other computations performed all use the lower precision. 

* **_gmres1_dq.m and gmres1.m_** are functions that run left-preconditioned GMRES using precisions single/double and double/quad.

* **_gmres2.m and gmres2_dq.m,_** are functions that perform Arnoldi iterations in GCRO-DR using precisions single/double and double/quad.

* **_getHarmVecs1.m and getHarmVecs2.m_** are functions in GCRO-DR library that determine harmonic Ritz vectors in GCRO-DR. The library and associated functions are available at https://www.sandia.gov/-mlparks/software/.

* **_gcrodr_dq.m, and gcrodr.m_** are functions that run GCRO-DR using precisions single/double and double/quad.

* **_rgmresir_test_randn.m_** is an example script for comparing GMRES-IR and RGMRES-IR (with 3 precisions) on random dense matrices with mode 3 having various condition numbers.

* **_rgmresir_test_ss.m_** is an example script for comparing GMRES-IR and RGMRES-IR (with 3 precisions) on matrices in SuiteSparse collection.


## Requirements
* The codes have been developed and tested with MATLAB 2020a.
* rgmresir_test_ss.m requires ssget MATLAB interface for testing the algorithm on matrices in SuiteSparse collection.
* The codes require the Advanpix Multiprecision Computing Toolbox for extended precision computations. 
A free trial of Advanpix is available for download from https://www.advanpix.com/.


