This folder contains the codes for numerical experiments A and B in the paper "Spectrum-Revealing Cholesky Factorization for Kernel Methods". 

files description:

abalone.txt and CCPP.txt are data sets. 

generate_distance_matrix.f90 and generate_RBF_kernel are used to create kernel matrices. 

partial_dgeqpf.f is modified dgeqpf.f. It can implement a partial QRCP.

srch.f90 implements SRCH with extra swaps.

random_projection_ch.f90 implements SRCH without extra swaps. 

svd_A_abalone.txt and svd_A_CCPP.txt store the singular values of kernel matrices of abalone and CCPP. Computing the singular values is time consuming. 

test_srch_A.f90 compares run time of SRCH and DPSTRF.

test_srch_B.f90 compares residual errors of SRCH and DPSTRF.

Implementation:
enter make, then ./job_test_srch_A.sh and ./job_test_srch_B.sh

the matlab codes to draw pictures are in folder graphs_AB_CCPP