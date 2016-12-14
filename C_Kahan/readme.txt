This folder contains the codes for numerical experiment C in the paper "Spectrum-Revealing Cholesky Factorization for Kernel Methods". 

files description:

abalone.txt and CCPP.txt are data sets. 

generate_distance_matrix.f90 and generate_RBF_kernel are used to create kernel matrices. 

partial_dgeqpf.f is modified dgeqpf.f. It can implement a partial QRCP.

srch.f90 implements SRCH with extra swaps.

random_projection_ch.f90 implements SRCH without extra swaps. 

svd_A_abalone.txt and svd_A_CCPP.txt store the singular values of kernel matrices of abalone and CCPP. Computing the singular values is time consuming. 

test_srch_C.f90 compares the results of SRCH and DPSTRF on Kahan matrix.

Implementation:
enter make, then ./test_srch_C