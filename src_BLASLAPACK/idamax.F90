       INTEGER FUNCTION idamax(N,DX,INCX)
 !     .. Scalar Arguments ..
       INTEGER incx,n
 !     .. Array Arguments ..
       DOUBLE PRECISION dx(*)
 !     .. Local Scalars ..
       DOUBLE PRECISION dmax
       INTEGER i,ix
 !     .. Intrinsic Functions ..
       INTRINSIC dabs
 
       idamax = 0
       IF (n.LT.1 .OR. incx.LE.0) RETURN
       idamax = 1
       IF (n.EQ.1) RETURN
       IF (incx.EQ.1) THEN
 
 !        code for increment equal to 1
 
          dmax = dabs(dx(1))
          DO i = 2,n
             IF (dabs(dx(i)).GT.dmax) THEN
                idamax = i
                dmax = dabs(dx(i))
             END IF
          END DO
       ELSE
 
 !        code for increment not equal to 1
 
          ix = 1
          dmax = dabs(dx(1))
          ix = ix + incx
          DO i = 2,n
             IF (dabs(dx(ix)).GT.dmax) THEN
                idamax = i
                dmax = dabs(dx(ix))
             END IF
             ix = ix + incx
          END DO
       END IF
       RETURN
       END