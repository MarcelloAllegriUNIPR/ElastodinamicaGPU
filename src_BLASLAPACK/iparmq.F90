       INTEGER FUNCTION iparmq( ISPEC, NAME, OPTS, N, ILO, IHI, LWORK )
 !     .. Scalar Arguments ..
       INTEGER            ihi, ilo, ispec, lwork, n
       CHARACTER          name*( * ), opts*( * )
 !     .. Parameters ..
       INTEGER            inmin, inwin, inibl, ishfts, iacc22
       parameter( inmin = 12, inwin = 13, inibl = 14, ishfts = 15, iacc22 = 16 )
       INTEGER            nmin, k22min, kacmin, nibble, knwswp
       parameter( nmin = 75, k22min = 14, kacmin = 14, nibble = 14, knwswp = 500 )
       REAL               two
       parameter( two = 2.0 )
 !     .. Local Scalars ..
       INTEGER            nh, ns
       INTEGER            i, ic, iz
       CHARACTER          subnam*6
 !     .. Intrinsic Functions ..
       INTRINSIC          log, max, mod, nint, real
 
 !     .. Executable Statements ..
       IF( ( ispec.EQ.ishfts ) .OR. ( ispec.EQ.inwin ) .OR. ( ispec.EQ.iacc22 ) ) THEN
 
 !        ==== Set the number simultaneous shifts ====
 
          nh = ihi - ilo + 1
          ns = 2
          IF( nh.GE.30 ) ns = 4
          IF( nh.GE.60 ) ns = 10
          IF( nh.GE.150 ) ns = max( 10, nh / nint( log( real( nh ) ) / log( two ) ) )
          IF( nh.GE.590 ) ns = 64
          IF( nh.GE.3000 ) ns = 128
          IF( nh.GE.6000 ) ns = 256
          ns = max( 2, ns-mod( ns, 2 ) )
       END IF
 
       IF( ispec.EQ.inmin ) THEN
 
 
 !        ===== Matrices of order smaller than NMIN get sent
 !        .     to xLAHQR, the classic double shift algorithm.
 !        .     This must be at least 11. ====
 
          iparmq = nmin
 
       ELSE IF( ispec.EQ.inibl ) THEN
 
 !        ==== INIBL: skip a multi-shift qr iteration and
 !        .    whenever aggressive early deflation finds
 !        .    at least (NIBBLE*(window size)/100) deflations. ====
 
          iparmq = nibble
 
       ELSE IF( ispec.EQ.ishfts ) THEN
 
 !        ==== NSHFTS: The number of simultaneous shifts =====
 
          iparmq = ns
 
       ELSE IF( ispec.EQ.inwin ) THEN
 
 !        ==== NW: deflation window size.  ====
 
          IF( nh.LE.knwswp ) THEN
             iparmq = ns
          ELSE
             iparmq = 3*ns / 2
          END IF
 
       ELSE IF( ispec.EQ.iacc22 ) THEN
 
 !        ==== IACC22: Whether to accumulate reflections
 !        .     before updating the far-from-diagonal elements
 !        .     and whether to use 2-by-2 block structure while
 !        .     doing it.  A small amount of work could be saved
 !        .     by making this choice dependent also upon the
 !        .     NH=IHI-ILO+1.
 
 
 !        Convert NAME to upper case if the first character is lower case.
 
          iparmq = 0
          subnam = name
          ic = ichar( subnam( 1: 1 ) )
          iz = ichar( 'Z' )
          IF( iz.EQ.90 .OR. iz.EQ.122 ) THEN
 
 !           ASCII character set
 
             IF( ic.GE.97 .AND. ic.LE.122 ) THEN
                subnam( 1: 1 ) = char( ic-32 )
                DO i = 2, 6
                   ic = ichar( subnam( i: i ) )
                   IF( ic.GE.97 .AND. ic.LE.122 ) subnam( i: i ) = char( ic-32 )
                END DO
             END IF
 
          ELSE IF( iz.EQ.233 .OR. iz.EQ.169 ) THEN
 
 !           EBCDIC character set
 
             IF( ( ic.GE.129 .AND. ic.LE.137 ) .OR. ( ic.GE.145 .AND. ic.LE.153 ) .OR. ( ic.GE.162 .AND. ic.LE.169 ) ) THEN
                subnam( 1: 1 ) = char( ic+64 )
                DO i = 2, 6
                   ic = ichar( subnam( i: i ) )
                   IF( ( ic.GE.129 .AND. ic.LE.137 ) .OR. ( ic.GE.145 .AND. ic.LE.153 ) .OR. ( ic.GE.162 .AND. ic.LE.169 ) )subnam( i:i ) = char( ic+64 )
                END DO
             END IF
 
          ELSE IF( iz.EQ.218 .OR. iz.EQ.250 ) THEN
 
 !           Prime machines:  ASCII+128
 
             IF( ic.GE.225 .AND. ic.LE.250 ) THEN
                subnam( 1: 1 ) = char( ic-32 )
                DO i = 2, 6
                   ic = ichar( subnam( i: i ) )
                   IF( ic.GE.225 .AND. ic.LE.250 ) subnam( i: i ) = char( ic-32 )
                END DO
             END IF
          END IF
 
          IF( subnam( 2:6 ).EQ.'GGHRD' .OR. subnam( 2:6 ).EQ.'GGHD3' ) THEN
             iparmq = 1
             IF( nh.GE.k22min ) iparmq = 2
          ELSE IF ( subnam( 4:6 ).EQ.'EXC' ) THEN
             IF( nh.GE.kacmin ) iparmq = 1
             IF( nh.GE.k22min ) iparmq = 2
          ELSE IF ( subnam( 2:6 ).EQ.'HSEQR' .OR. subnam( 2:5 ).EQ.'LAQR' ) THEN
             IF( ns.GE.kacmin ) iparmq = 1
             IF( ns.GE.k22min ) iparmq = 2
          END IF
 
       ELSE
 !        ===== invalid value of ispec =====
          iparmq = -1
 
       END IF
 
       END