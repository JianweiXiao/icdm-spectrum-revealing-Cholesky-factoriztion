This folder contains the fortran codes to compute Gaussian process on CCPP data matrix. This folder contains the codes for numerical experiment E in the paper "Spectrum-Revealing Cholesky Factorization for Kernel Methods". 
We compare the performance on SRCH and DPSTRF. For the details of the prediction formulas used. Check the paper "Stable and Efficient Gaussian Process Calculations".

files description:

CCPP_test_data_6568_4.txt:	test data matrix
CCPP_train_data_3000_4.txt:	train data matrix
CCPP_test_label_6568_1.txt: labels of test data matrix
CCPP_train_label_3000_1.txt: labels of train data matrix

generate_distance_matrix.f90, generate_distance_matrix1_matrix2.f90, generate_RBF_kernel.f90, generate_RBF_kernel_matrix1_matrix2.f90 are routines to compute covariance matrix. 

partial_dgeqpf.f: modified version of partial_dgeqpf, can implement a partial QRCP

srch.f90: spectrum revealing cholesky factorization 

test_srch_E.f90: main program, compare the Gaussian process of SRCH and DPSTRF on the CCPP data

job_test_srch_E.sh: script on different tolerance 

Implementation:
enter make, then ./job_test_srch_E.sh

the matlab codes to draw pictures are in folder graphs_E_GP