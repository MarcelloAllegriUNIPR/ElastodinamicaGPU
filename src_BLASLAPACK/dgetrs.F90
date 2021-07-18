       SUBROUTINE dgetrs( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
 !     .. Scalar Arguments ..
       CHARACTER          TRANS
       INTEGER            INFO, LDA, LDB, N, NRHS
 !     .. Array Arguments ..
       INTEGER            IPIV( * )
       DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
 !     .. Parameters ..
       DOUBLE PRECISION   ONE
       parameter( one = 1.0d+0 )
 !     .. Local Scalars ..
       LOGICAL            NOTRAN
 !     .. External Functions ..
       LOGICAL            LSAME
       EXTERNAL           lsame
 !     .. External Subroutines ..
       EXTERNAL           dlaswp, dtrsm, xerbla
 !     .. Intrinsic Functions ..
       INTRINSIC          max
 
 !     .. Executable Statements ..
 
 !     Test the input parameters.
 
       info = 0
       notran = lsame( trans, 'N' )
       IF( .NOT.notran .AND. .NOT.lsame( trans, 'T' ) .AND. .NOT. lsame( trans, 'C' ) ) THEN
          info = -1
       ELSE IF( n.LT.0 ) THEN
          info = -2
       ELSE IF( nrhs.LT.0 ) THEN
          info = -3
       ELSE IF( lda.LT.max( 1, n ) ) THEN
          info = -5
       ELSE IF( ldb.LT.max( 1, n ) ) THEN
          info = -8
       END IF
       IF( info.NE.0 ) THEN
          CALL xerbla( 'DGETRS', -info )
          RETURN
       END IF
 
 !     Quick return if possible
 
       IF( n.EQ.0 .OR. nrhs.EQ.0 ) RETURN
 
       IF( notran ) THEN
 
 !        Solve A * X = B.
 
 !        Apply row interchanges to the right hand sides.
 
          CALL dlaswp( nrhs, b, ldb, 1, n, ipiv, 1 )
 
 !        Solve L*X = B, overwriting B with X.
 
          CALL dtrsm( 'Left', 'Lower', 'No transpose', 'Unit', n, nrhs, one, a, lda, b, ldb )
 
 !        Solve U*X = B, overwriting B with X.
 
          CALL dtrsm( 'Left', 'Upper', 'No transpose', 'Non-unit', n, nrhs, one, a, lda, b, ldb )
       ELSE
 
 !        Solve A**T * X = B.
 
 !        Solve U**T *X = B, overwriting B with X.
 
          CALL dtrsm( 'Left', 'Upper', 'Transpose', 'Non-unit', n, nrhs, one, a, lda, b, ldb )
 
 !        Solve L**T *X = B, overwriting B with X.
 
          CALL dtrsm( 'Left', 'Lower', 'Transpose', 'Unit', n, nrhs, one, a, lda, b, ldb )
 
 !        Apply row interchanges to the solution vectors.
 
          CALL dlaswp( nrhs, b, ldb, 1, n, ipiv, -1 )
       END IF
 
       RETURN
 
       END