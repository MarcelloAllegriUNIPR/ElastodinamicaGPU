       INTEGER FUNCTION ilaenv( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
 !     .. Scalar Arguments ..
       CHARACTER*( * )    name, opts
       INTEGER            ispec, n1, n2, n3, n4
 !     .. Local Scalars ..
       INTEGER            i, ic, iz, nb, nbmin, nx
       LOGICAL            cname, sname, twostage
       CHARACTER          c1*1, c2*2, c4*2, c3*3, subnam*16
 !     .. Intrinsic Functions ..
       INTRINSIC          char, ichar, int, min, real
 !     .. External Functions ..
       INTEGER            ieeeck, iparmq, iparam2stage
       EXTERNAL           ieeeck, iparmq, iparam2stage
 
 !     .. Executable Statements ..
 
       GO TO ( 10, 10, 10, 80, 90, 100, 110, 120,130, 140, 150, 160, 160, 160, 160, 160)ispec
 
 !     Invalid value for ISPEC
 
       ilaenv = -1
       RETURN
 
    10 CONTINUE
 
 !     Convert NAME to upper case if the first character is lower case.
 
       ilaenv = 1
       subnam = name
       ic = ichar( subnam( 1: 1 ) )
       iz = ichar( 'Z' )
       IF( iz.EQ.90 .OR. iz.EQ.122 ) THEN
 
 !        ASCII character set
 
          IF( ic.GE.97 .AND. ic.LE.122 ) THEN
             subnam( 1: 1 ) = char( ic-32 )
             DO 20 i = 2, 6
                ic = ichar( subnam( i: i ) )
                IF( ic.GE.97 .AND. ic.LE.122 )subnam( i: i ) = char( ic-32 )
    20       CONTINUE
          END IF
 
       ELSE IF( iz.EQ.233 .OR. iz.EQ.169 ) THEN
 
 !        EBCDIC character set
 
          IF( ( ic.GE.129 .AND. ic.LE.137 ) .OR. ( ic.GE.145 .AND. ic.LE.153 ) .OR. ( ic.GE.162 .AND. ic.LE.169 ) ) THEN
             subnam( 1: 1 ) = char( ic+64 )
             DO 30 i = 2, 6
                ic = ichar( subnam( i: i ) )
                IF( ( ic.GE.129 .AND. ic.LE.137 ) .OR. ( ic.GE.145 .AND. ic.LE.153 ) .OR. ( ic.GE.162 .AND. ic.LE.169 ) )subnam( i:i ) = char( ic+64 )
    30       CONTINUE
          END IF
 
       ELSE IF( iz.EQ.218 .OR. iz.EQ.250 ) THEN
 
 !        Prime machines:  ASCII+128
 
          IF( ic.GE.225 .AND. ic.LE.250 ) THEN
             subnam( 1: 1 ) = char( ic-32 )
             DO 40 i = 2, 6
                ic = ichar( subnam( i: i ) )
                IF( ic.GE.225 .AND. ic.LE.250 ) subnam( i: i ) = char( ic-32 )
    40       CONTINUE
          END IF
       END IF
 
       c1 = subnam( 1: 1 )
       sname = c1.EQ.'S' .OR. c1.EQ.'D'
       cname = c1.EQ.'C' .OR. c1.EQ.'Z'
       IF( .NOT.( cname .OR. sname ) ) RETURN
       c2 = subnam( 2: 3 )
       c3 = subnam( 4: 6 )
       c4 = c3( 2: 3 )
       twostage = len( subnam ).GE.11.AND. subnam( 11: 11 ).EQ.'2'
 
       GO TO ( 50, 60, 70 )ispec
 
    50 CONTINUE
 
 !     ISPEC = 1:  block size
 
 !     In these examples, separate code is provided for setting NB for
 !     real and complex.  We assume that NB will take the same value in
 !     single or double precision.
 
       nb = 1
 
       IF( subnam(2:6).EQ.'LAORH' ) THEN
 
 !        This is for *LAORHR_GETRFNP routine
 
          IF( sname ) THEN
              nb = 32
          ELSE
              nb = 32
          END IF
       ELSE IF( c2.EQ.'GE' ) THEN
          IF( c3.EQ.'TRF' ) THEN
             IF( sname ) THEN
                nb = 64
             ELSE
                nb = 64
             END IF
          ELSE IF( c3.EQ.'QRF' .OR. c3.EQ.'RQF' .OR. c3.EQ.'LQF' .OR.c3.EQ.'QLF' ) THEN
             IF( sname ) THEN
                nb = 32
             ELSE
                nb = 32
             END IF
          ELSE IF( c3.EQ.'QR ') THEN
             IF( n3 .EQ. 1) THEN
                IF( sname ) THEN
 !     M*N
                   IF ((n1*n2.LE.131072).OR.(n1.LE.8192)) THEN
                      nb = n1
                   ELSE
                      nb = 32768/n2
                   END IF
                ELSE
                   IF ((n1*n2.LE.131072).OR.(n1.LE.8192)) THEN
                      nb = n1
                   ELSE
                      nb = 32768/n2
                   END IF
                END IF
             ELSE
                IF( sname ) THEN
                   nb = 1
                ELSE
                   nb = 1
                END IF
             END IF
          ELSE IF( c3.EQ.'LQ ') THEN
             IF( n3 .EQ. 2) THEN
                IF( sname ) THEN
 !     M*N
                   IF ((n1*n2.LE.131072).OR.(n1.LE.8192)) THEN
                      nb = n1
                   ELSE
                      nb = 32768/n2
                   END IF
                ELSE
                   IF ((n1*n2.LE.131072).OR.(n1.LE.8192)) THEN
                      nb = n1
                   ELSE
                      nb = 32768/n2
                   END IF
                END IF
             ELSE
                IF( sname ) THEN
                   nb = 1
                ELSE
                   nb = 1
                END IF
             END IF
          ELSE IF( c3.EQ.'HRD' ) THEN
             IF( sname ) THEN
                nb = 32
             ELSE
                nb = 32
             END IF
          ELSE IF( c3.EQ.'BRD' ) THEN
             IF( sname ) THEN
                nb = 32
             ELSE
                nb = 32
             END IF
          ELSE IF( c3.EQ.'TRI' ) THEN
             IF( sname ) THEN
                nb = 64
             ELSE
                nb = 64
             END IF
          END IF
       ELSE IF( c2.EQ.'PO' ) THEN
          IF( c3.EQ.'TRF' ) THEN
             IF( sname ) THEN
                nb = 64
             ELSE
                nb = 64
             END IF
          END IF
       ELSE IF( c2.EQ.'SY' ) THEN
          IF( c3.EQ.'TRF' ) THEN
             IF( sname ) THEN
                IF( twostage ) THEN
                   nb = 192
                ELSE
                   nb = 64
                END IF
             ELSE
                IF( twostage ) THEN
                   nb = 192
                ELSE
                   nb = 64
                END IF
             END IF
          ELSE IF( sname .AND. c3.EQ.'TRD' ) THEN
             nb = 32
          ELSE IF( sname .AND. c3.EQ.'GST' ) THEN
             nb = 64
          END IF
       ELSE IF( cname .AND. c2.EQ.'HE' ) THEN
          IF( c3.EQ.'TRF' ) THEN
             IF( twostage ) THEN
                nb = 192
             ELSE
                nb = 64
             END IF
          ELSE IF( c3.EQ.'TRD' ) THEN
             nb = 32
          ELSE IF( c3.EQ.'GST' ) THEN
             nb = 64
          END IF
       ELSE IF( sname .AND. c2.EQ.'OR' ) THEN
          IF( c3( 1: 1 ).EQ.'G' ) THEN
             IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ. 'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' )THEN
                nb = 32
             END IF
          ELSE IF( c3( 1: 1 ).EQ.'M' ) THEN
             IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ. 'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' )THEN
                nb = 32
             END IF
          END IF
       ELSE IF( cname .AND. c2.EQ.'UN' ) THEN
          IF( c3( 1: 1 ).EQ.'G' ) THEN
             IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ. 'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' ) THEN
                nb = 32
             END IF
          ELSE IF( c3( 1: 1 ).EQ.'M' ) THEN
             IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ. 'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' ) THEN
                nb = 32
             END IF
          END IF
       ELSE IF( c2.EQ.'GB' ) THEN
          IF( c3.EQ.'TRF' ) THEN
             IF( sname ) THEN
                IF( n4.LE.64 ) THEN
                   nb = 1
                ELSE
                   nb = 32
                END IF
             ELSE
                IF( n4.LE.64 ) THEN
                   nb = 1
                ELSE
                   nb = 32
                END IF
             END IF
          END IF
       ELSE IF( c2.EQ.'PB' ) THEN
          IF( c3.EQ.'TRF' ) THEN
             IF( sname ) THEN
                IF( n2.LE.64 ) THEN
                   nb = 1
                ELSE
                   nb = 32
                END IF
             ELSE
                IF( n2.LE.64 ) THEN
                   nb = 1
                ELSE
                   nb = 32
                END IF
             END IF
          END IF
       ELSE IF( c2.EQ.'TR' ) THEN
          IF( c3.EQ.'TRI' ) THEN
             IF( sname ) THEN
                nb = 64
             ELSE
                nb = 64
             END IF
          ELSE IF ( c3.EQ.'EVC' ) THEN
             IF( sname ) THEN
                nb = 64
             ELSE
                nb = 64
             END IF
          END IF
       ELSE IF( c2.EQ.'LA' ) THEN
          IF( c3.EQ.'UUM' ) THEN
             IF( sname ) THEN
                nb = 64
             ELSE
                nb = 64
             END IF
          END IF
       ELSE IF( sname .AND. c2.EQ.'ST' ) THEN
          IF( c3.EQ.'EBZ' ) THEN
             nb = 1
          END IF
       ELSE IF( c2.EQ.'GG' ) THEN
          nb = 32
          IF( c3.EQ.'HD3' ) THEN
             IF( sname ) THEN
                nb = 32
             ELSE
                nb = 32
             END IF
          END IF
       END IF
       ilaenv = nb
       RETURN
 
    60 CONTINUE
 
 !     ISPEC = 2:  minimum block size
 
       nbmin = 2
       IF( c2.EQ.'GE' ) THEN
          IF( c3.EQ.'QRF' .OR. c3.EQ.'RQF' .OR. c3.EQ.'LQF' .OR. c3.EQ. 'QLF' ) THEN
             IF( sname ) THEN
                nbmin = 2
             ELSE
                nbmin = 2
             END IF
          ELSE IF( c3.EQ.'HRD' ) THEN
             IF( sname ) THEN
                nbmin = 2
             ELSE
                nbmin = 2
             END IF
          ELSE IF( c3.EQ.'BRD' ) THEN
             IF( sname ) THEN
                nbmin = 2
             ELSE
                nbmin = 2
             END IF
          ELSE IF( c3.EQ.'TRI' ) THEN
             IF( sname ) THEN
                nbmin = 2
             ELSE
                nbmin = 2
             END IF
          END IF
       ELSE IF( c2.EQ.'SY' ) THEN
          IF( c3.EQ.'TRF' ) THEN
             IF( sname ) THEN
                nbmin = 8
             ELSE
                nbmin = 8
             END IF
          ELSE IF( sname .AND. c3.EQ.'TRD' ) THEN
             nbmin = 2
          END IF
       ELSE IF( cname .AND. c2.EQ.'HE' ) THEN
          IF( c3.EQ.'TRD' ) THEN
             nbmin = 2
          END IF
       ELSE IF( sname .AND. c2.EQ.'OR' ) THEN
          IF( c3( 1: 1 ).EQ.'G' ) THEN
             IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ. 'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' ) THEN
                nbmin = 2
             END IF
          ELSE IF( c3( 1: 1 ).EQ.'M' ) THEN
             IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ. 'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' )THEN
                nbmin = 2
             END IF
          END IF
       ELSE IF( cname .AND. c2.EQ.'UN' ) THEN
          IF( c3( 1: 1 ).EQ.'G' ) THEN
             IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ. 'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' ) THEN
                nbmin = 2
             END IF
          ELSE IF( c3( 1: 1 ).EQ.'M' ) THEN
             IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ. 'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' )THEN
                nbmin = 2
             END IF
          END IF
       ELSE IF( c2.EQ.'GG' ) THEN
          nbmin = 2
          IF( c3.EQ.'HD3' ) THEN
             nbmin = 2
          END IF
       END IF
       ilaenv = nbmin
       RETURN
 
    70 CONTINUE
 
 !     ISPEC = 3:  crossover point
 
       nx = 0
       IF( c2.EQ.'GE' ) THEN
          IF( c3.EQ.'QRF' .OR. c3.EQ.'RQF' .OR. c3.EQ.'LQF' .OR. c3.EQ. 'QLF' ) THEN
             IF( sname ) THEN
                nx = 128
             ELSE
                nx = 128
             END IF
          ELSE IF( c3.EQ.'HRD' ) THEN
             IF( sname ) THEN
                nx = 128
             ELSE
                nx = 128
             END IF
          ELSE IF( c3.EQ.'BRD' ) THEN
             IF( sname ) THEN
                nx = 128
             ELSE
                nx = 128
             END IF
          END IF
       ELSE IF( c2.EQ.'SY' ) THEN
          IF( sname .AND. c3.EQ.'TRD' ) THEN
             nx = 32
          END IF
       ELSE IF( cname .AND. c2.EQ.'HE' ) THEN
          IF( c3.EQ.'TRD' ) THEN
             nx = 32
          END IF
       ELSE IF( sname .AND. c2.EQ.'OR' ) THEN
          IF( c3( 1: 1 ).EQ.'G' ) THEN
             IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ. 'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' ) THEN
                nx = 128
             END IF
          END IF
       ELSE IF( cname .AND. c2.EQ.'UN' ) THEN
          IF( c3( 1: 1 ).EQ.'G' ) THEN
             IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ. 'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' ) THEN
                nx = 128
             END IF
          END IF
       ELSE IF( c2.EQ.'GG' ) THEN
          nx = 128
          IF( c3.EQ.'HD3' ) THEN
             nx = 128
          END IF
       END IF
       ilaenv = nx
       RETURN
 
    80 CONTINUE
 
 !     ISPEC = 4:  number of shifts (used by xHSEQR)
 
       ilaenv = 6
       RETURN
 
    90 CONTINUE
 
 !     ISPEC = 5:  minimum column dimension (not used)
 
       ilaenv = 2
       RETURN
 
   100 CONTINUE
 
 !     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
 
       ilaenv = int( real( min( n1, n2 ) )*1.6e0 )
       RETURN
 
   110 CONTINUE
 
 !     ISPEC = 7:  number of processors (not used)
 
       ilaenv = 1
       RETURN
 
   120 CONTINUE
 
 !     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
 
       ilaenv = 50
       RETURN
 
   130 CONTINUE
 
 !     ISPEC = 9:  maximum size of the subproblems at the bottom of the
 !                 computation tree in the divide-and-conquer algorithm
 !                 (used by xGELSD and xGESDD)
 
       ilaenv = 25
       RETURN
 
   140 CONTINUE
 
 !     ISPEC = 10: ieee NaN arithmetic can be trusted not to trap
 
 !     ILAENV = 0
       ilaenv = 1
       IF( ilaenv.EQ.1 ) THEN
          ilaenv = ieeeck( 1, 0.0, 1.0 )
       END IF
       RETURN
 
   150 CONTINUE
 
 !     ISPEC = 11: infinity arithmetic can be trusted not to trap
 
 !     ILAENV = 0
       ilaenv = 1
       IF( ilaenv.EQ.1 ) THEN
          ilaenv = ieeeck( 0, 0.0, 1.0 )
       END IF
       RETURN
 
   160 CONTINUE
 
 !     12 <= ISPEC <= 16: xHSEQR or related subroutines.
 
       ilaenv = iparmq( ispec, name, opts, n1, n2, n3, n4 )
       RETURN
 
       END