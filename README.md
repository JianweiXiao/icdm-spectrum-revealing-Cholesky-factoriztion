# icdm-spectrum-revealing-Cholesky-factoriztion

These four folders contain the fortran or matlab codes for numerical experiments in the paper 
Spectrum-Revealing Cholesky Factorization for Kernel Methods. Jianwei Xiao and Ming Gu. IEEE International Conference on Data Mining (ICDM), Barcelona, Spain, 2016.

*********************************************************************************************

Folder “AB_CCPP” contains the fortran code for numerical experiments A and B (CCPP matrix). 

Implementation:
enter make, then ./job_test_srch_A.sh and ./job_test_srch_B.sh

the matlab codes to draw pictures are in folder graphs_AB_CCPP. SRCH without swap will already get a good low rank approximation.

*********************************************************************************************

Folder “C_Kahan” contains the fortran code for numerical experiment C (Kahan matrix). 

Implementation:
enter make, then ./test_srch_C

It will print out the singular value approximation ratio. The first column is for dpstrf, the second column is for SRCH without swap, the third column is for SRCH.

*********************************************************************************************

Folder “D_CSI” contains the matlab code for numerical experiment D (comparison with CSI). 

Implementation:
run my_comparison_digit.m, you will get the figure 4 and figure 5. To get figure 6, modify the parameter m to [50,100,150,200] and obtain the run time. To get figure 6, run running_time_csi_srch.m

*********************************************************************************************

Folder “E_GP” contains the fortran code for numerical experiment E (Gaussian process). 

Implementation:
enter make, then ./job_test_srch_E.sh

the matlab codes to draw pictures are in folder graphs_E_GP

*********************************************************************************************

Citing this work

We ask those who benefit from this work to cite the paper 
Spectrum-Revealing Cholesky Factorization for Kernel Methods. Jianwei Xiao and Ming Gu. IEEE International Conference on Data Mining (ICDM), Barcelona, Spain, 2016.