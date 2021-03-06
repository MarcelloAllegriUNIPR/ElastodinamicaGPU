       SUBROUTINE dgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
 !     .. Scalar Arguments ..
       DOUBLE PRECISION ALPHA,BETA
       INTEGER K,LDA,LDB,LDC,M,N
       CHARACTER TRANSA,TRANSB
 !     .. Array Arguments ..
       DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
 !     .. External Functions ..
       LOGICAL LSAME
       EXTERNAL lsame
 !     .. External Subroutines ..
       EXTERNAL xerbla
 !     .. Intrinsic Functions ..
       INTRINSIC max
 !     .. Local Scalars ..
       DOUBLE PRECISION TEMP
       INTEGER I,INFO,J,L,NCOLA,NROWA,NROWB
       LOGICAL NOTA,NOTB
 !     .. Parameters ..
       DOUBLE PRECISION ONE,ZERO
       parameter(one=1.0d+0,zero=0.0d+0)
 
 !     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
 !     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
 !     and  columns of  A  and the  number of  rows  of  B  respectively.
 
       nota = lsame(transa,'N')
       notb = lsame(transb,'N')
       IF (nota) THEN
           nrowa = m
           ncola = k
       ELSE
           nrowa = k
           ncola = m
       END IF
       IF (notb) THEN
           nrowb = k
       ELSE
           nrowb = n
       END IF
 
 !     Test the input parameters.
 
       info = 0
       IF ((.NOT.nota) .AND. (.NOT.lsame(transa,'C')) .AND. (.NOT.lsame(transa,'T'))) THEN
           info = 1
       ELSE IF ((.NOT.notb) .AND. (.NOT.lsame(transb,'C')) .AND. (.NOT.lsame(transb,'T'))) THEN
           info = 2
       ELSE IF (m.LT.0) THEN
           info = 3
       ELSE IF (n.LT.0) THEN
           info = 4
       ELSE IF (k.LT.0) THEN
           info = 5
       ELSE IF (lda.LT.max(1,nrowa)) THEN
           info = 8
       ELSE IF (ldb.LT.max(1,nrowb)) THEN
           info = 10
       ELSE IF (ldc.LT.max(1,m)) THEN
           info = 13
       END IF
       IF (info.NE.0) THEN
           CALL xerbla('DGEMM ',info)
           RETURN
       END IF
 
 !     Quick return if possible.
 
       IF ((m.EQ.0) .OR. (n.EQ.0) .OR.(((alpha.EQ.zero).OR. (k.EQ.0)).AND. (beta.EQ.one))) RETURN
 
 !     And if  alpha.eq.zero.
 
       IF (alpha.EQ.zero) THEN
           IF (beta.EQ.zero) THEN
               DO 20 j = 1,n
                   DO 10 i = 1,m
                       c(i,j) = zero
    10             CONTINUE
    20         CONTINUE
           ELSE
               DO 40 j = 1,n
                   DO 30 i = 1,m
                       c(i,j) = beta*c(i,j)
    30             CONTINUE
    40         CONTINUE
           END IF
           RETURN
       END IF
 
 !     Start the operations.
 
       IF (notb) THEN
           IF (nota) THEN
 
 !           Form  C := alpha*A*B + beta*C.
 
               DO 90 j = 1,n
                   IF (beta.EQ.zero) THEN
                       DO 50 i = 1,m
                           c(i,j) = zero
    50                 CONTINUE
                   ELSE IF (beta.NE.one) THEN
                       DO 60 i = 1,m
                           c(i,j) = beta*c(i,j)
    60                 CONTINUE
                   END IF
                   DO 80 l = 1,k
                       temp = alpha*b(l,j)
                       DO 70 i = 1,m
                           c(i,j) = c(i,j) + temp*a(i,l)
    70                 CONTINUE
    80             CONTINUE
    90         CONTINUE
           ELSE
 
 !           Form  C := alpha*A**T*B + beta*C
 
               DO 120 j = 1,n
                   DO 110 i = 1,m
                       temp = zero
                       DO 100 l = 1,k
                           temp = temp + a(l,i)*b(l,j)
   100                 CONTINUE
                       IF (beta.EQ.zero) THEN
                           c(i,j) = alpha*temp
                       ELSE
                           c(i,j) = alpha*temp + beta*c(i,j)
                       END IF
   110             CONTINUE
   120         CONTINUE
           END IF
       ELSE
           IF (nota) THEN
 
 !           Form  C := alpha*A*B**T + beta*C
 
               DO 170 j = 1,n
                   IF (beta.EQ.zero) THEN
                       DO 130 i = 1,m
                           c(i,j) = zero
   130                 CONTINUE
                   ELSE IF (beta.NE.one) THEN
                       DO 140 i = 1,m
                           c(i,j) = beta*c(i,j)
   140                 CONTINUE
                   END IF
                   DO 160 l = 1,k
                       temp = alpha*b(j,l)
                       DO 150 i = 1,m
                           c(i,j) = c(i,j) + temp*a(i,l)
   150                 CONTINUE
   160             CONTINUE
   170         CONTINUE
           ELSE
 
 !           Form  C := alpha*A**T*B**T + beta*C
 
               DO 200 j = 1,n
                   DO 190 i = 1,m
                       temp = zero
                       DO 180 l = 1,k
                           temp = temp + a(l,i)*b(j,l)
   180                 CONTINUE
                       IF (beta.EQ.zero) THEN
                           c(i,j) = alpha*temp
                       ELSE
                           c(i,j) = alpha*temp + beta*c(i,j)
                       END IF
   190             CONTINUE
   200         CONTINUE
           END IF
       END IF
 
       RETURN
 
       END