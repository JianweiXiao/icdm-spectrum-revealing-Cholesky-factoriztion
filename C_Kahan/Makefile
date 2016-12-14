CC = gfortran
CFLAGS = 
LIBS = -O3 -lblas -llapack
TARGETS	= test_srch_C

all: 	$(TARGETS)

# executable files

test_srch_C: test_srch_C.o generate_distance_matrix.o generate_RBF_kernel.o random_projection_ch.o srch.o dlapmr.o partial_dgeqpf.o
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

test_srch_C.o: test_srch_C.f90
		$(CC) -c $(CFLAGS) $^ $(LIBS)	

# lapack subroutines
dlapmr.o: dlapmr.f
		$(CC) -c $(CFLAGS) $^ $(LIBS)	


clean :
	rm -f *.o
	rm -r $(TARGETS)
	rm -f *.out

