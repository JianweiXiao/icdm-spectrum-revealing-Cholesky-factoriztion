! MAIN PROGRAM FOR NUMERICAL EXPERIMENT A IN THE PAPER "Spectrum-Revealing Cholesky Factorization for Kernel Methods"

! COMPARE DPSTRF, SRCH AND RANDOM_PROJECTION_CH

! DPSTRF IS LAPACK SUBROUTINE
! RANDOM_PROJECTION_CH IS SRCH WITHOUT ADDING THE SWAP STRATEGY, SRCH NO SWAP
! SRCH 

PROGRAM TEST_SRCH_A

   IMPLICIT NONE 

   CHARACTER(10) ::    R_STRING, TOL_STRING
   INTEGER, PARAMETER :: DIM = 9568
   INTEGER, PARAMETER :: PARAMETER_NP = 30
   INTEGER, PARAMETER :: PARAMETER_NB = 20
   INTEGER, PARAMETER :: PARAMETER_R = 20
   INTEGER, PARAMETER :: COUNTER = 4
   REAL*8, PARAMETER :: SIGMA = 1.0
   REAL*8, PARAMETER :: G = 2
   REAL*8 :: TOL
   REAL*8, PARAMETER :: REGULARIZAR = 1e-9
   INTEGER, PARAMETER :: ITERATION = 1

   INTEGER :: N = DIM
   INTEGER :: NP = PARAMETER_NP
   INTEGER :: NB = PARAMETER_NB
   INTEGER :: R
   REAL*8, DIMENSION(:,:), ALLOCATABLE :: A, B, C, W, W1, W2, W3, A1, A2, A3, L1, L2, L3, U, VT, AKBK, AK
   INTEGER, DIMENSION(:,:), ALLOCATABLE :: PC1, PC2, PC3, PR1, PR2, PR3
   REAL*8, DIMENSION(:), ALLOCATABLE :: AKBKS1, AKBKS2, AKBKS3, AS, AKS1, AKS2, AKS3
   INTEGER :: ISEED(4)
   INTEGER :: I, J, INFO
   INTEGER, DIMENSION(:), ALLOCATABLE :: PIV 
   INTEGER, DIMENSION(:), ALLOCATABLE :: WORK
   REAL*8 :: DLANGE, A1NORM
   REAL*8, DIMENSION(:,:), ALLOCATABLE :: RESIDUAL1, RESIDUAL2, RESIDUAL3
   REAL*8, DIMENSION(:,:), ALLOCATABLE :: ORIGINAL
   REAL*8, PARAMETER :: ONE = 1.0_8
   REAL*8, PARAMETER :: ZERO = 1.0_8
   REAL*8, DIMENSION(:), ALLOCATABLE :: LOG_SVD
   INTEGER :: M, LWORK
   REAL*8, ALLOCATABLE, DIMENSION(:) :: W_WORK
   INTEGER, ALLOCATABLE, DIMENSION(:) :: IWORK, IFAIL
   REAL*8, ALLOCATABLE, DIMENSION(:,:) :: Z, TEMP_A
   REAL*8 :: RHO, PHI, RHO1, PHI1, KESAI, KESAI1
   REAL*8, ALLOCATABLE, DIMENSION(:,:) :: S, K
   INTEGER :: RANK
   REAL*8 :: START, FINISH

   CALL GETARG(1, R_STRING)
   READ(R_STRING, *) R
   WRITE (*,*) "APPROXIMATE RANK IS", R
   CALL GETARG(2, TOL_STRING)
   READ(TOL_STRING, *) TOL
   WRITE (*,*) "TOLERANCE IS", TOL

   ISEED(1) = 40
   ISEED(2) = 829
   ISEED(3) = 90
   ISEED(4) = 65

   ! RANDOM TEST MATRIX
   IF (COUNTER == 1) THEN
      WRITE (*,*) 'RANDOM MATRIX'
      WRITE (*,*) 'N =', N
      WRITE (*,*) 'R =', R
      WRITE (*,*) 'NB =', NB
      WRITE (*,*) 'NP =', NP
      ALLOCATE(A(N,N), A1(N,N), A2(N,N), A3(N,N), ORIGINAL(N,N))
      CALL DLARNV(3, ISEED, N*N, A)
      DO I = 1,N
         DO J = I+1,N
            A(I,J) = 0.0_8
         END DO
      END DO
      A = MATMUL(A, TRANSPOSE(A))
      A = 0.5 * (A + TRANSPOSE(A))
   ELSE
   END IF

   ! KAHAN MATRIX
   IF (COUNTER == 2) THEN
      WRITE (*,*) 'KAHAN MATRIX'
      WRITE (*,*) 'N =', N
      WRITE (*,*) 'R =', R
      WRITE (*,*) 'NB =', NB
      WRITE (*,*) 'NP =', NP
      ALLOCATE(A(N,N), A1(N,N), A2(N,N), A3(N,N), ORIGINAL(N,N))
      ALLOCATE(S(N,N), K(N,N))
      S = 0.0_8
      K = 0.0_8
      RHO = 0.9999
      PHI1 = 0.285
      KESAI1 = SQRT(1.0-PHI1*PHI1)
      PHI = PHI1*RHO
      KESAI = KESAI1*RHO
      DO I = 1,N
         S(I,I) = KESAI ** (I-1)
      END DO
      DO I = 1,N
         K(I,I) = 1.0_8
         IF (I < N) THEN
            DO J = I+1,N
               K(I,J) = -PHI
            END DO
         ELSE
         END IF
      END DO
      A = MATMUL(S,K)
      A = MATMUL(TRANSPOSE(A),A)
      A = 0.5 * (TRANSPOSE(A) + A)
      DEALLOCATE(S, K)
   ELSE
   END IF


   ! REAL DATA: ABALONE
   IF (COUNTER == 3) THEN
      WRITE (*,*) 'ABALONE'
      N = 4177
      R = PARAMETER_R
      NB = PARAMETER_NB
      NP = PARAMETER_NP
      WRITE (*,*) 'N =', N
      WRITE (*,*) 'R =', R
      WRITE (*,*) 'NB =', NB
      WRITE (*,*) 'NP =', NP
      ALLOCATE(A(N,N), A1(N,N), A2(N,N), A3(N,N), ORIGINAL(N,N))
      ALLOCATE(B(4177,8), C(4177,4177))

      OPEN(UNIT = 1, FILE = "abalone.txt", ACTION = "READ")
      DO I = 1,4177
         READ(1,*) (B(I,J), J = 1,8)
      END DO

      CALL GENERATE_DISTANCE_MATRIX(B, 4177, 8, C)
      CALL GENERATE_RBF_KERNEL(C, 4177, SIGMA, A)
      A = 0.5 * (A + TRANSPOSE(A))

      ! ADD REGULARIZATION 
      DO I = 1,N
         A(I,I) = A(I,I) + REGULARIZAR
      END DO
      DEALLOCATE(B, C)
   ELSE
   END IF

   ! REAL DATA: CCPP
   IF (COUNTER == 4) THEN
      WRITE (*,*) 'CCPP'
      N = 9568
      NB = PARAMETER_NB
      NP = PARAMETER_NP
      WRITE (*,*) 'N =', N
      WRITE (*,*) 'R =', R
      WRITE (*,*) 'NB =', NB
      WRITE (*,*) 'NP =', NP
      ALLOCATE(A(N,N), A1(N,N), A2(N,N), A3(N,N), ORIGINAL(N,N))
      ALLOCATE(B(9568,4), C(9568,9568))

      OPEN(UNIT = 1, FILE = "CCPP.txt", ACTION = "READ")
      DO I = 1,9568
         READ(1,*) (B(I,J), J = 1,4)
      END DO

      CALL GENERATE_DISTANCE_MATRIX(B, 9568, 4, C)
      CALL GENERATE_RBF_KERNEL(C, 9568, SIGMA, A)
      A = 0.5 * (A + TRANSPOSE(A))

      ! ADD REGULARIZATION 
      DO I = 1,N
         A(I,I) = A(I,I) + REGULARIZAR
      END DO
      DEALLOCATE(B, C)
   ELSE
   END IF

   ALLOCATE(W(NP,N), W1(NP,N), W2(NP,N), W3(NP,N), L1(N,N), L2(N,N), L3(N,N), PC1(1,N), &
      PC2(1,N), PC3(1,N), PR1(1,N), PR2(1,N), PR3(1,N), WORK(10*N))
   ISEED(1) = 54
   ISEED(2) = 98
   ISEED(2) = 67
   ISEED(4) = 197
   CALL DLARNV(3, ISEED, NP*N, W) 

   ORIGINAL = A
   A1 = A
   A2 = A
   A3 = A
   W1 = W
   W2 = W
   W3 = W

   ALLOCATE(AKBK(N,R), U(N,N), VT(N,N), AKBKS1(N), AKBKS2(N), AKBKS3(N), AS(N), AKS1(N), AKS2(N), &
      AKS3(N), AK(R,R))
   AKBKS1 = 0.0_8
   AKBKS2 = 0.0_8
   AKBKS3 = 0.0_8
   AKS1 = 0.0_8
   AKS2 = 0.0_8
   AKS3 = 0.0_8
   AK = 0.0_8
   AS = 0.0_8
   U = 0.0_8
   VT = 0.0_8

   IF (COUNTER == 3) THEN
      ! ONLY RUN ONE TO STORE THE SINGULAR VALUES OF A IN TXT FILE
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !   CALL DGESVD('N', 'N', N, N, ORIGINAL, N, AS, U, N, VT, N, WORK, 5*N, INFO)
      !   WRITE (*,*) 'COMPUTE SINGULAR VALUES OF A FINISHED'
      !   OPEN(UNIT = 1, FILE = "SVD_A_ABALONE.TXT")
      !   DO I = 1,N
      !      WRITE (1,*) AS(I) 
      !   END DO
      !   CLOSE(1)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      OPEN(UNIT = 1, FILE = "svd_A_abalone.txt", ACTION = "READ")
      DO I = 1,N
         READ(1,*) AS(I)
      END DO
      CLOSE(1)

   ELSE IF (COUNTER == 4) THEN
      ! ONLY RUN ONE TO STORE THE SINGULAR VALUES OF A IN TXT FILE
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !CHOOSE REGULARIZATION = 1E-9
      !CALL DGESVD('N', 'N', N, N, ORIGINAL, N, AS, U, N, VT, N, WORK, 5*N, INFO)
      !WRITE (*,*) 'COMPUTE SINGULAR VALUES OF A FINISHED'
      !OPEN(UNIT = 2, FILE = "SVD_A_CCPP.TXT")
      !DO I = 1,N
      !   WRITE (2,*) AS(I)  
      !END DO
      !CLOSE(2)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      OPEN(UNIT = 2, FILE = "svd_A_CCPP.txt", ACTION = "READ")
      DO I = 1,N
         READ(2,*) AS(I)
      END DO
      CLOSE(2)
   ELSE
      CALL DGESVD('N', 'N', N, N, ORIGINAL, N, AS, U, N, VT, N, WORK, 5*N, INFO)
      WRITE (*,*) 'COMPUTE SINGULAR VALUES OF A FINISHED'
   END IF

   ! WRITE (*,*) '***************************************'
   ALLOCATE(PIV(N))
   A1 = A
   CALL CPU_TIME(START)
   CALL DPSTRF('L', N, A1, N, PIV, RANK, TOL, WORK, INFO)
   CALL CPU_TIME(FINISH)
   WRITE (*,*) 'DPSTRF'
   WRITE (*,*) 'TIME =', (FINISH - START)/ITERATION
   WRITE (*,*) 'RANK OF DPSTRF', RANK
   L1 = 0.0_8
   DO I = 1,N
      DO J = 1,MIN(I,R)
         L1(I,J) = A1(I,J)
      END DO
   END DO
   PC1(1,:) = PIV
   PR1(1,:) = PIV
   DEALLOCATE(PIV)
   ! WRITE (*,*) '***************************************'

   ! WRITE (*,*) '***************************************'
   ALLOCATE(PIV(N))
   ! WRITE (*,*) 'RUN SRCH NO SWAP' 
   W2 = W
   A2 = A
   CALL CPU_TIME(START)
   CALL RANDOM_PROJECTION_CH(A2, N, W2, NB, NP, R, PIV)
   CALL CPU_TIME(FINISH)
   WRITE (*,*) 'SRCH NO SWAP'
   WRITE (*,*) 'TIME =', (FINISH - START)/ITERATION
   L2 = 0.0_8
   DO I = 1,N
      DO J = 1,MIN(I,R)
         L2(I,J) = A2(I,J)
      END DO
   END DO
   PC2(1,:) = PIV
   PR2(1,:) = PIV
   DEALLOCATE(PIV)
   ! WRITE (*,*) '***************************************'

   ! WRITE (*,*) '***************************************'
   ALLOCATE(PIV(N))
   ! WRITE (*,*) 'RUN SRCH'
   W3 = W
   A3 = A
   CALL CPU_TIME(START)
   CALL SRCH(A3, N, W3, NB, NP, R, PIV, G)
   CALL CPU_TIME(FINISH)
   ! WRITE (*,*) 'SRCH'
   ! WRITE (*,*) 'TIME =', (FINISH - START)/ITERATION
   L3 = 0.0_8
   DO I = 1,N
      DO J = 1,MIN(I,R)
         L3(I,J) = A3(I,J)
      END DO
   END DO
   PC3(1,:) = PIV
   PR3(1,:) = PIV
   DEALLOCATE(PIV)
   ! WRITE (*,*) '***************************************'
   WRITE (*,*)

   ! COMPUTE SINGULAR VALUES OF AKBK1, AKBK2, AKBK3
   ! AKBK = L1(1:N,1:R)
   ! CALL DGESVD('N', 'N', N, R, AKBK, N, AKBKS1, U, N, VT, N, WORK, 5*N, INFO)
   ! AKBK = L2(1:N,1:R)
   ! CALL DGESVD('N', 'N', N, R, AKBK, N, AKBKS2, U, N, VT, N, WORK, 5*N, INFO)
   ! AKBK = L3(1:N,1:R)
   ! CALL DGESVD('N', 'N', N, R, AKBK, N, AKBKS3, U, N, VT, N, WORK, 5*N, INFO)

   ! ! COMPUTE SINGULAR VALUES OF AK1, AK2, AK3
   ! AK = L1(1:R,1:R)
   ! CALL DGESVD('N', 'N', R, R, AK, R, AKS1, U, N, VT, N, WORK, 5*N, INFO)
   ! AK = L2(1:R,1:R)
   ! CALL DGESVD('N', 'N', R, R, AK, R, AKS2, U, N, VT, N, WORK, 5*N, INFO)
   ! AK = L3(1:R,1:R)
   ! CALL DGESVD('N', 'N', R, R, AK, R, AKS3, U, N, VT, N, WORK, 5*N, INFO)

   ! WRITE (*,*) 'RESIDUAL ERROR, (SIGMA_J(A)-SIGMA_J(AKBK)^2)/SIGMA_J(A)'
   ! WRITE (*,'(20G12.4)') 'INDEX', 'DPSTRF ', 'SRCH NO SWAP'
   ! DO I = 1,10
   !    WRITE (*,'(20G12.4)') I, (AS(I)-AKBKS1(I)**2)/AS(I), (AS(I)-AKBKS2(I)**2)/AS(I)
   ! END DO
   ! WRITE (*,*)

   ! WRITE (*,*) 'SIGMA_J(AKBK)^2/SIGMA_J(A)'
   ! WRITE (*,'(20G12.4)') '     INDEX ', 'DPSTRF ', 'SRCH NO SWAP ', 'SRCH'
   ! DO I = 1,10
   !    WRITE (*,'(20G12.4)') I, (AKBKS1(I)**2)/AS(I), (AKBKS2(I)**2)/AS(I), (AKBKS3(I)**2)/AS(I)
   ! END DO
   ! WRITE (*,*)

   ! CHECK CORRECTNESS
   ! ALLOCATE(RESIDUAL1(N,N), RESIDUAL2(N,N), RESIDUAL3(N,N))
   ! RESIDUAL1 = 0.0_8
   ! RESIDUAL2 = 0.0_8
   ! RESIDUAL3 = 0.0_8
   ! CALL DGEMM('N', 'T', N, N, N, ONE, L1, N, L1, N, ZERO, RESIDUAL1, N)
   ! CALL DGEMM('N', 'T', N, N, N, ONE, L2, N, L2, N, ZERO, RESIDUAL2, N)
   ! CALL DGEMM('N', 'T', N, N, N, ONE, L3, N, L3, N, ZERO, RESIDUAL3, N)
   ! RESIDUAL1 = A(PR1(1,:), PC1(1,:)) - RESIDUAL1
   ! RESIDUAL2 = A(PR2(1,:), PC2(1,:)) - RESIDUAL2
   ! RESIDUAL3 = A(PR3(1,:), PC3(1,:)) - RESIDUAL3

   ! WRITE (*,*) 'DPSTRF, SHOULD BE CLOSE TO ZERO'
   ! WRITE (*,*) MAXVAL(RESIDUAL1(1:R,1:N))
   ! WRITE (*,*) MAXVAL(-RESIDUAL1(1:R,1:N))
   ! WRITE (*,*) 'SRCH WITHOUT SWAP STRATEGY, SHOULD BE CLOSE TO ZERO'
   ! WRITE (*,*) MAXVAL(RESIDUAL2(1:R,1:N))
   ! WRITE (*,*) MAXVAL(-RESIDUAL2(1:R,1:N))
   ! WRITE (*,*) 'SRCH, SHOULD BE ZERO'
   ! WRITE (*,*) MAXVAL(RESIDUAL3(1:R,1:N))
   ! WRITE (*,*) MAXVAL(-RESIDUAL3(1:R,1:N))

   ! DEALLOCATE(RESIDUAL1, RESIDUAL2, RESIDUAL3, ORIGINAL)
   
   DEALLOCATE(W, W1, W2, W3, L1, L2, L3, PC1, PC2, PC3, PR1, PR2, PR3, WORK)
   DEALLOCATE(AKBK, U, VT, AKBKS1, AKBKS2, AKBKS3, AS, AK)
   DEALLOCATE(AKS1, AKS2, AKS3)
END PROGRAM TEST_SRCH_A
