CC = gfortran
CFLAGS = 
LIBS = -O3 -lblas -llapack
TARGETS	= test_srch_A test_srch_B

all: 	$(TARGETS)

# executable files

test_srch_A: test_srch_A.o generate_distance_matrix.o generate_RBF_kernel.o random_projection_ch.o srch.o partial_dgeqpf.o dlapmr.o
		$(CC) $^ -o $@ $(LIBS)

test_srch_B: test_srch_B.o generate_distance_matrix.o generate_RBF_kernel.o random_projection_ch.o srch.o partial_dgeqpf.o dlapmr.o
		$(CC) $^ -o $@ $(LIBS)

# objective files

generate_distance_matrix.o: generate_distance_matrix.f90
		$(CC) -c $(CFLAGS) $^ $(LIBS)

generate_RBF_kernel.o: generate_RBF_kernel.f90
		$(CC) -c $(CFLAGS) $^ $(LIBS)				

random_projection_ch.o: random_projection_ch.f90
		$(CC) -c $(CFLAGS) $^ $(LIBS)

srch.o: srch.f90
		$(CC) -c $(CFLAGS) $^ $(LIBS)

partial_dgeqpf.o: partial_dgeqpf.f
		$(CC) -c $(CFLAGS) $^ $(LIBS)

# test objective

test_srch_A.o: test_srch_A.f90
		$(CC) -c $(CFLAGS) $^ $(LIBS)	

test_srch_B.o: test_srch_B.f90
		$(CC) -c $(CFLAGS) $^ $(LIBS)	

# lapack routines

dlapmr.o: dlapmr.f
		$(CC) -c $(CFLAGS) $^ $(LIBS)	

clean :
	rm -f *.o
	rm -r $(TARGETS)
	rm -f *.out

