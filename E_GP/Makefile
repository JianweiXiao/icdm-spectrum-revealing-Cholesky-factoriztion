CC = gfortran
CFLAGS = 
LIBS = -O3 -lblas -llapack
TARGETS	= test_srch_E

all: 	$(TARGETS)

# executable files
test_srch_E: test_srch_E.o random_projection_ch.o generate_distance_matrix.o generate_RBF_kernel.o generate_distance_matrix1_matrix2.o generate_RBF_kernel_matrix1_matrix2.o partial_dgeqpf.o dlapmr.o
		$(CC) $^ -o $@ $(LIBS)

# objective files
random_projection_ch.o: random_projection_ch.f90
		$(CC) -c $(CFLAGS) $^ $(LIBS)

generate_distance_matrix.o: generate_distance_matrix.f90
		$(CC) -c $(CFLAGS) $^ $(LIBS)

generate_RBF_kernel.o: generate_RBF_kernel.f90
		$(CC) -c $(CFLAGS) $^ $(LIBS)

generate_distance_matrix1_matrix2.o: generate_distance_matrix1_matrix2.f90
		$(CC) -c $(CFLAGS) $^ $(LIBS)

generate_RBF_kernel_matrix1_matrix2.o: generate_RBF_kernel_matrix1_matrix2.f90
		$(CC) -c $(CFLAGS) $^ $(LIBS)

partial_dgeqpf.o: partial_dgeqpf.f
		$(CC) -c $(CFLAGS) $^ $(LIBS)

# test objective
test_srch_E.o: test_srch_E.f90	
		$(CC) -c $(CFLAGS) $^ $(LIBS)

# lapack routines

dlapmr.o: dlapmr.f
		$(CC) -c $(CFLAGS) $^ $(LIBS)	

clean :
	rm -f *.o
	rm -r $(TARGETS)
	rm -f *.out

