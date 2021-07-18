       SUBROUTINE dgetrf( M, N, A, LDA, IPIV, INFO )
!     .. Scalar Arguments ..
       INTEGER            INFO, LDA, M, N
!     .. Array Arguments ..
       INTEGER            IPIV( * )
       DOUBLE PRECISION   A( LDA, * )
!     .. Parameters ..
       DOUBLE PRECISION   ONE
       parameter( one = 1.0d+0 )
!     .. Local Scalars ..
       INTEGER            I, IINFO, J, JB, NB
!     .. External Subroutines ..
       EXTERNAL           dgemm, dgetrf2, dlaswp, dtrsm, xerbla
!     .. External Functions ..
       INTEGER            ILAENV
       EXTERNAL           ilaenv
!     .. Intrinsic Functions ..
       INTRINSIC          max, min
       
!     .. Executable Statements ..
!     Test the input parameters.

       info = 0
       IF( m.LT.0 ) THEN
          info = -1
       ELSE IF( n.LT.0 ) THEN
          info = -2
       ELSE IF( lda.LT.max( 1, m ) ) THEN
          info = -4
       END IF
       IF( info.NE.0 ) THEN
          CALL xerbla( 'DGETRF', -info )
          RETURN
       END IF
!     Quick return if possible
       IF( m.EQ.0 .OR. n.EQ.0 ) RETURN
!     Determine the block size for this environment.
       nb = ilaenv( 1, 'DGETRF', ' ', m, n, -1, -1 )
       IF( nb.LE.1 .OR. nb.GE.min( m, n ) ) THEN
!        Use unblocked code.
          CALL dgetrf2( m, n, a, lda, ipiv, info )
       ELSE
!        Use blocked code.
          DO 20 j = 1, min( m, n ), nb
             jb = min( min( m, n )-j+1, nb )
!           Factor diagonal and subdiagonal blocks and test for exact
!           singularity.
             CALL dgetrf2( m-j+1, jb, a( j, j ), lda, ipiv( j ), iinfo )
!           Adjust INFO and the pivot indices.
             IF( info.EQ.0 .AND. iinfo.GT.0 )info = iinfo + j - 1
             DO 10 i = j, min( m, j+jb-1 )
                ipiv( i ) = j - 1 + ipiv( i )
    10       CONTINUE
!           Apply interchanges to columns 1:J-1.
             CALL dlaswp( j-1, a, lda, j, j+jb-1, ipiv, 1 )
             IF( j+jb.LE.n ) THEN
!              Apply interchanges to columns J+JB:N.
                CALL dlaswp( n-j-jb+1, a( 1, j+jb ), lda, j, j+jb-1,ipiv, 1 )
!              Compute block row of U.
                CALL dtrsm( 'Left', 'Lower', 'No transpose', 'Unit', jb,n-j-jb+1, one, a( j, j ), lda, a( j, j+jb ),lda )
                IF( j+jb.LE.m ) THEN
!                 Update trailing submatrix.
                   CALL dgemm( 'No transpose', 'No transpose', m-j-jb+1,n-j-jb+1, jb, -one, a( j+jb, j ), lda,a( j, j+jb ), lda, one, a( j+jb, j+jb ),lda )
                END IF
             END IF
    20    CONTINUE
       END IF
       RETURN
       END