       RECURSIVE SUBROUTINE dgetrf2( M, N, A, LDA, IPIV, INFO )
!     .. Scalar Arguments ..
       INTEGER            info, lda, m, n
!    .. Array Arguments ..
       INTEGER            ipiv( * )
       DOUBLE PRECISION   a( lda, * )
!     .. Parameters ..
       DOUBLE PRECISION   one, zero
       parameter( one = 1.0d+0, zero = 0.0d+0 )
!     .. Local Scalars ..
       DOUBLE PRECISION   sfmin, temp
       INTEGER            i, iinfo, n1, n2
!     .. External Functions ..
       DOUBLE PRECISION   dlamch
       INTEGER            idamax
       EXTERNAL           dlamch, idamax
!     .. External Subroutines ..
       EXTERNAL           dgemm, dscal, dlaswp, dtrsm, xerbla
!     .. Intrinsic Functions ..
       INTRINSIC          max, min
       
!     .. Executable Statements ..
!     Test the input parameters

       info = 0
       IF( m.LT.0 ) THEN
          info = -1
       ELSE IF( n.LT.0 ) THEN
          info = -2
       ELSE IF( lda.LT.max( 1, m ) ) THEN
          info = -4
       END IF
       IF( info.NE.0 ) THEN
          CALL xerbla( 'DGETRF2', -info )
          RETURN
       END IF
 
!     Quick return if possible

       IF( m.EQ.0 .OR. n.EQ.0 )RETURN
  
       IF ( m.EQ.1 ) THEN
!        Use unblocked code for one row case
!        Just need to handle IPIV and INFO

          ipiv( 1 ) = 1
          IF ( a(1,1).EQ.zero )info = 1

       ELSE IF( n.EQ.1 ) THEN
!        Use unblocked code for one column case
!        Compute machine safe minimum
          sfmin = dlamch('S')
!        Find pivot and test for singularity
          i = idamax( m, a( 1, 1 ), 1 )
          ipiv( 1 ) = i
          IF( a( i, 1 ).NE.zero ) THEN
!           Apply the interchange
             IF( i.NE.1 ) THEN
                temp = a( 1, 1 )
                a( 1, 1 ) = a( i, 1 )
                a( i, 1 ) = temp
             END IF
!           Compute elements 2:M of the column
             IF( abs(a( 1, 1 )) .GE. sfmin ) THEN
                CALL dscal( m-1, one / a( 1, 1 ), a( 2, 1 ), 1 )
             ELSE
                DO 10 i = 1, m-1
                   a( 1+i, 1 ) = a( 1+i, 1 ) / a( 1, 1 )
    10          CONTINUE
             END IF

          ELSE
             info = 1
          END IF

       ELSE

 !        Use recursive code

          n1 = min( m, n ) / 2
          n2 = n-n1
 
 !               [ A11 ]
 !        Factor [ --- ]
 !               [ A21 ]
 
          CALL dgetrf2( m, n1, a, lda, ipiv, iinfo )
  
          IF ( info.EQ.0 .AND. iinfo.GT.0 )info = iinfo
 !                              [ A12 ]
 !        Apply interchanges to [ --- ]
 !                              [ A22 ]
          CALL dlaswp( n2, a( 1, n1+1 ), lda, 1, n1, ipiv, 1 )
 !        Solve A12
          CALL dtrsm( 'L', 'L', 'N', 'U', n1, n2, one, a, lda,a( 1, n1+1 ), lda )
 !        Update A22
          CALL dgemm( 'N', 'N', m-n1, n2, n1, -one, a( n1+1, 1 ), lda,a( 1, n1+1 ), lda, one, a( n1+1, n1+1 ), lda )
 !        Factor A22
          CALL dgetrf2( m-n1, n2, a( n1+1, n1+1 ), lda, ipiv( n1+1 ),iinfo )
 !        Adjust INFO and the pivot indices
          IF ( info.EQ.0 .AND. iinfo.GT.0 )info = iinfo + n1
          DO 20 i = n1+1, min( m, n )
             ipiv( i ) = ipiv( i ) + n1
    20    CONTINUE
 !        Apply interchanges to A21
          CALL dlaswp( n1, a( 1, 1 ), lda, n1+1, min( m, n), ipiv, 1 )
 
       END IF
       RETURN
       END