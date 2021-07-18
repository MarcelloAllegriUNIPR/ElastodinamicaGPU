       SUBROUTINE dlaswp( N, A, LDA, K1, K2, IPIV, INCX )
 !     .. Scalar Arguments ..
       INTEGER            INCX, K1, K2, LDA, N
 !     .. Array Arguments ..
       INTEGER            IPIV( * )
       DOUBLE PRECISION   A( LDA, * )
 !     .. Local Scalars ..
       INTEGER            I, I1, I2, INC, IP, IX, IX0, J, K, N32
       DOUBLE PRECISION   TEMP
 
 !     .. Executable Statements ..
 !     Interchange row I with row IPIV(K1+(I-K1)*abs(INCX)) for each of rows
 !     K1 through K2.
 
       IF( incx.GT.0 ) THEN
          ix0 = k1
          i1 = k1
          i2 = k2
          inc = 1
       ELSE IF( incx.LT.0 ) THEN
          ix0 = k1 + ( k1-k2 )*incx
          i1 = k2
          i2 = k1
          inc = -1
       ELSE
          RETURN
       END IF
 
       n32 = ( n / 32 )*32
       IF( n32.NE.0 ) THEN
          DO 30 j = 1, n32, 32
             ix = ix0
             DO 20 i = i1, i2, inc
                ip = ipiv( ix )
                IF( ip.NE.i ) THEN
                   DO 10 k = j, j + 31
                      temp = a( i, k )
                      a( i, k ) = a( ip, k )
                      a( ip, k ) = temp
    10             CONTINUE
                END IF
                ix = ix + incx
    20       CONTINUE
    30    CONTINUE
       END IF
       IF( n32.NE.n ) THEN
          n32 = n32 + 1
          ix = ix0
          DO 50 i = i1, i2, inc
             ip = ipiv( ix )
             IF( ip.NE.i ) THEN
                DO 40 k = n32, n
                   temp = a( i, k )
                   a( i, k ) = a( ip, k )
                   a( ip, k ) = temp
    40          CONTINUE
             END IF
             ix = ix + incx
    50    CONTINUE
       END IF
 
       RETURN
 
       END