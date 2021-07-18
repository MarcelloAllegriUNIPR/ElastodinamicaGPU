       SUBROUTINE xerbla( SRNAME, INFO )
 !     .. Scalar Arguments ..
       CHARACTER*(*)      SRNAME
       INTEGER            INFO
 !     .. Intrinsic Functions ..
       INTRINSIC          len_trim
 
 !     .. Executable Statements ..
 
      WRITE( *, fmt = 9999 )srname( 1:len_trim( srname ) ), info
 
      stop
 
9999  FORMAT( ' ** On entry to ', a, ' parameter number ', i2, ' had ','an illegal value' )
 
      END