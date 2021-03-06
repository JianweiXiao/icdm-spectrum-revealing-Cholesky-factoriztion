! COMPUTE DISTANCE MATRIX USING ONE DATA MATRIX

SUBROUTINE GENERATE_DISTANCE_MATRIX(DATAPOINTS, M, N, DISTMAT)

! DISTMAT = GENERATE_DISTANCE_MATRIX(DATAPOINTS)
!   
! THE ROWS OF DATAPOINTS ARE THE OBSERVATIONS X_I
! RETURNS DISTMAT SO THE (IJ)TH ENTRY IS NORM(X_I-X_J)^2
!
! FIRST WHITENS DATA: CENTERS AND MAKES STD DEV OF POINTS 1 
! BEFORE CALCULATIONS

INTEGER, INTENT(IN) :: M, N
REAL*8, INTENT(INOUT), DIMENSION(M,N) :: DATAPOINTS
REAL*8, INTENT(OUT), DIMENSION(M,M) :: DISTMAT

REAL*8, DIMENSION(N) :: MEAN ! MEAN OF EVERY COLUMN OF DATAPOINTS
REAL*8, DIMENSION(N) :: SUM_SQUARE
INTEGER :: I,J,K

DISTMAT = 0.0_8
SUM_SQUARE = 0.0_8
MEAN = SUM(DATAPOINTS, 1) / M 

! SUBTRACT MEAN ON EVERY COLUMN OF DATAPOINTS
DO I = 1,N
   DATAPOINTS(:,I) = DATAPOINTS(:,I) - MEAN(I)
END DO

DO J = 1,N
   DO I = 1,M
      SUM_SQUARE(J) = SUM_SQUARE(J) + DATAPOINTS(I,J)**2
   END DO
END DO

SUM_SQUARE = SQRT(SUM_SQUARE/(M-1))
DO I = 1,N
   DATAPOINTS(:,I) = DATAPOINTS(:,I) / SUM_SQUARE(I)
END DO

! NOW BUILD THE DISTANCE MATRIX. DISTMAT(I,J) = ||XI-XJ||_2^2
DO I = 1,M
   DO J = 1,I-1
      DO K = 1,N
         DISTMAT(I,J) = DISTMAT(I,J) + (DATAPOINTS(I,K)-DATAPOINTS(J,K))**2
      END DO
      DISTMAT(J,I) = DISTMAT(I,J) 
   END DO
END DO

END SUBROUTINE GENERATE_DISTANCE_MATRIX
