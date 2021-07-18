       INTEGER FUNCTION ieeeck( ISPEC, ZERO, ONE )
 !     .. Scalar Arguments ..
       INTEGER            ispec
       REAL               one, zero
 !     .. Local Scalars ..
       REAL               nan1, nan2, nan3, nan4, nan5, nan6, neginf, negzro, newzro, posinf
 
 !     .. Executable Statements ..
       ieeeck = 1
 
       posinf = one / zero
       IF( posinf.LE.one ) THEN
          ieeeck = 0
          RETURN
       END IF
 
       neginf = -one / zero
       IF( neginf.GE.zero ) THEN
          ieeeck = 0
          RETURN
       END IF
 
       negzro = one / ( neginf+one )
       IF( negzro.NE.zero ) THEN
          ieeeck = 0
          RETURN
       END IF
 
       neginf = one / negzro
       IF( neginf.GE.zero ) THEN
          ieeeck = 0
          RETURN
       END IF
 
       newzro = negzro + zero
       IF( newzro.NE.zero ) THEN
          ieeeck = 0
          RETURN
       END IF
 
       posinf = one / newzro
       IF( posinf.LE.one ) THEN
          ieeeck = 0
          RETURN
       END IF
 
       neginf = neginf*posinf
       IF( neginf.GE.zero ) THEN
          ieeeck = 0
          RETURN
       END IF
 
       posinf = posinf*posinf
       IF( posinf.LE.one ) THEN
          ieeeck = 0
          RETURN
       END IF
 
 !     Return if we were only asked to check infinity arithmetic
 
       IF( ispec.EQ.0 ) RETURN
 
       nan1 = posinf + neginf
 
       nan2 = posinf / neginf
 
       nan3 = posinf / posinf
 
       nan4 = posinf*zero
 
       nan5 = neginf*negzro
 
       nan6 = nan5*zero
 
       IF( nan1.EQ.nan1 ) THEN
          ieeeck = 0
          RETURN
       END IF
 
       IF( nan2.EQ.nan2 ) THEN
          ieeeck = 0
          RETURN
       END IF
 
       IF( nan3.EQ.nan3 ) THEN
          ieeeck = 0
          RETURN
       END IF
 
       IF( nan4.EQ.nan4 ) THEN
          ieeeck = 0
          RETURN
       END IF
 
       IF( nan5.EQ.nan5 ) THEN
          ieeeck = 0
          RETURN
       END IF
 
       IF( nan6.EQ.nan6 ) THEN
          ieeeck = 0
          RETURN
       END IF
 
       RETURN
       END