       DOUBLE PRECISION FUNCTION ddot(N,DX,INCX,DY,INCY)
 !     .. Scalar Arguments ..
       INTEGER incx,incy,n
 !     .. Array Arguments ..
       DOUBLE PRECISION dx(*),dy(*)
 !     .. Local Scalars ..
       DOUBLE PRECISION dtemp
       INTEGER i,ix,iy,m,mp1
 !     .. Intrinsic Functions ..
       INTRINSIC mod
       ddot = 0.0d0
       dtemp = 0.0d0
       IF (n.LE.0) RETURN
       IF (incx.EQ.1 .AND. incy.EQ.1) THEN
 
 !        code for both increments equal to 1
 
 !        clean-up loop
 
          m = mod(n,5)
          IF (m.NE.0) THEN
             DO i = 1,m
                dtemp = dtemp + dx(i)*dy(i)
             END DO
             IF (n.LT.5) THEN
                ddot=dtemp
             RETURN
             END IF
          END IF
          mp1 = m + 1
          DO i = mp1,n,5
           dtemp = dtemp + dx(i)*dy(i) + dx(i+1)*dy(i+1) + dx(i+2)*dy(i+2) + dx(i+3)*dy(i+3) + dx(i+4)*dy(i+4)
          END DO
       ELSE
 
 !        code for unequal increments or equal increments
 !          not equal to 1
 
          ix = 1
          iy = 1
          IF (incx.LT.0) ix = (-n+1)*incx + 1
          IF (incy.LT.0) iy = (-n+1)*incy + 1
          DO i = 1,n
             dtemp = dtemp + dx(ix)*dy(iy)
             ix = ix + incx
             iy = iy + incy
          END DO
       END IF
       ddot = dtemp
       RETURN
       END