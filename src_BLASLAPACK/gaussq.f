	SUBROUTINE GAUSSQ(KIND,N,ALPHA,BETA,KPTS,ENDPTS,B,T,W)
C
C
C       Questa routine determina nodi e pesi di formule
C       -----------------------------------------------
C       di quadratura Gaussiane classiche
C       ---------------------------------
C
C       La routine e' scritta in FORTRAN IV, doppia precisione
C
C       Finalita' della routine
C       -----------------------
C
C       La routine calcola i nodi T(J) e i pesi W(J) di formule
C       di quadratura Gaussiane classiche, incluse quelle di tipo
C       Radau e Lobatto, che possono essere utilizzate per ap-
C       prossimare integrali del tipo
C
C                    I     W(x)F(x)dx
C                     (a,b)
C
C       dove W(x) puo' assumere una delle sei forme descritte piu'
C       avanti.
C
C       Metodo utilizzato
C       -----------------
C       I coefficienti della relazione di ricorrenza a tre termini
C       caratterizzante i polinomi ortogonali in questione vengono
C       utilizzati per costruire una matrice tridiagonale simmetri-
C       ca i cui autovalori, determinati con il metodo QL implicito
C       con shift, coincidono con i nodi della formula Gaussiana
C       cercata. La prima componente di ciascun autovettore orto-
C       normalizzato fornira' invece il peso della formula associa-
C       to al nodo (autovalore) corrispondente.
C
C       Bibliografia
C       ------------
C       G.H.Golub, J.H.Welsch, " Calculation of Gaussian quadrature
C             rules", Math.Comp., v.23, 1969, pp.221-230.
C       G.H.Golub, "Some modified matrix eigenvalue problems", SIAM
C             Review, v.15, 1973, pp.318-344.
C       G.Monegato, Calcolo Numerico, Levrotto & Bella, Torino, 1985.
C
C       Lista parametri
C       ---------------
C       In ingresso:
C
C       KIND - intero da 1 a 6 che indica il tipo di formula scelto:
C              KIND=1: formula di Legendre, W(x)=1 in (-1,1)
C              KIND=2: formula di Chebyshev di prima specie,
C                      W(x)=1/SQRT(1-x*x) in (-1,1)
C              KIND=3: formula di Chebyshev di seconda specie,
C                      W(x)=SQRT(1-x*x) in (-1,1)
C              KIND=4: formula di Hermite, W(x)=EXP(-x*x) in (-inf,+inf)
C              KIND=5: formula di Jacobi, W(x)=(1-x)**ALPHA*(1+x)**BETA
C                      in (-1,1), ALPHA e BETA > -1
C              KIND=6: formula di Laguerre generalizzata, W(x)=EXP(-x)
C                      in (0,inf)
C       N    - numero complessivo di nodi della formula scelta
C       ALPHA- parametro reale da inserire nei casi KIND=5 e KIND=6
C       BETA - parametro reale da inserire solo nel caso KIND=5
C       KPTS - intero normalmente posto =0, tranne quando uno dei
C              due estremi dell'intervallo (-1,1) (caso di Radau, KPTS=1)
C              oppure entrambi (caso di Lobatto, KPTS=2) sono nodi
C       ENDPTS-vettore reale di lunghezza 2 contenente gli eventuali
C              estremi -1., oppure 1., oppure -1.,1. scelti come nodi
C              (casi KPTS=1 e KPTS=2)
C       B    - vettore reale di servizio di lunghezza N
C
C       In uscita:
C
C       T    - vettore reale di lunghezza N contenente i nodi della
C              formula scelta
C       W    - vettore reale di lunghezza N contenente i pesi della
C              formula scelta
C
C       Sottoprogrammi richiamati
C       -------------------------
C              SOLVE
C              CLASS
C              IMTQL2
C              GAMMA
C              DSQRT, DBLE, FLOAT, DABS, DSIGN (funzioni FORTRAN di libreria)
C
	DOUBLE PRECISION B,T,W,ENDPTS,MUZERO,T1,GAM,SOLVE,DSQRT,ALPHA,BETA
	DIMENSION B(N),T(N),W(N),ENDPTS(2)
	CALL CLASS(KIND,N,ALPHA,BETA,B,T,MUZERO)
	IF(KPTS.EQ.0) GOTO 100
	IF(KPTS.EQ.2) GOTO 50
	T(N)=SOLVE(ENDPTS(1),N,T,B)*B(N-1)**2+ENDPTS(1)
	GOTO 100
50      GAM=SOLVE(ENDPTS(1),N,T,B)
	T1=((ENDPTS(1)-ENDPTS(2))/(SOLVE(ENDPTS(2),N,T,B)-GAM))
	B(N-1)=DSQRT(T1)
	T(N)=ENDPTS(1)+GAM*T1
100     W(1)=1.D0
	DO 105 I=2,N
105     W(I)=0.D0
	CALL IMTQL2(N,T,B,W,IERR)
	DO 110 I=1,N
110     W(I)=MUZERO*W(I)*W(I)
	RETURN
	END
C
	DOUBLE PRECISION FUNCTION SOLVE(SHIFT,N,A,B)
	REAL*8 SHIFT,A,B,ALPHA
	DIMENSION A(N),B(N)
	ALPHA=A(1)-SHIFT
	NM1=N-1
	DO 10 I=2,NM1
10      ALPHA=A(I)-SHIFT-B(I-1)**2/ALPHA
	SOLVE=1.D0/ALPHA
	RETURN
	END
C
	SUBROUTINE CLASS(KIND,N,ALPHA,BETA,B,A,MUZERO)
C
C       Scopo di questa subroutine e' la determinazione dei coefficienti
C       della relazione di ricorrenza a tre termini, memorizzati nei
C       vettori A e B, definita dai polinomi ortogonali scelti
C
	DOUBLE PRECISION A,B,MUZERO,ALPHA,BETA,ABI,A2B2,PI,GAMMA,DSQRT,AB
	DIMENSION A(N),B(N)
	DATA PI/3.141592653589793D0/
	NM1=N-1
	GO TO(10,20,30,40,50,60),KIND
10      MUZERO=2.D0
	DO 11 I=1,NM1
	A(I)=0.D0
	ABI=I
11      B(I)=ABI/DSQRT(4.D0*ABI*ABI-1.D0)
	A(N)=0.D0
	GO TO 70
20      MUZERO=PI
	DO 21 I=1,NM1
	A(I)=0.D0
21      B(I)=.5D0
	B(1)=DSQRT(0.5D0)
	A(N)=0.D0
	GOTO 70
30      MUZERO=PI/2.D0
	DO 31 I=1,NM1
	A(I)=0.D0
31      B(I)=0.5D0
	A(N)=0.D0
	GOTO 70
40      MUZERO=DSQRT(PI)
	DO 41 I=1,NM1
	A(I)=0.D0
41      B(I)=DSQRT(DBLE(FLOAT(I))/2.D0)
	A(N)=0.D0
	GOTO 70
50      AB=ALPHA+BETA
	ABI=2.D0+AB
	MUZERO=2.D0**(AB+1.D0)*GAMMA(ALPHA+1.D0)*GAMMA(BETA+1.D0)/GAMMA
     *  (ABI)
	A(1)=(BETA-ALPHA)/ABI
	B(1)=DSQRT(4.D0*(1.D0+ALPHA)*(1.D0+BETA)/((ABI+1.D0)*ABI*ABI))
	A2B2=BETA*BETA-ALPHA*ALPHA
	DO 51 I=2,NM1
	ABI=2.D0*DBLE(FLOAT(I))+AB
	A(I)=A2B2/((ABI-2.D0)*ABI)
51      B(I)=DSQRT(4.D0*DBLE(FLOAT(I))*(DBLE(FLOAT(I))+ALPHA)*(DBLE(
     *  FLOAT(I))+BETA)*(DBLE(FLOAT(I))+AB)/((ABI*ABI-1.D0)*ABI*ABI))
	ABI=2.D0*DBLE(FLOAT(N))+AB
	A(N)=A2B2/((ABI-2.D0)*ABI)
	GOTO 70
60      MUZERO=GAMMA(ALPHA+1.D0)
	DO 61 I=1,NM1
	A(I)=2.D0*DBLE(FLOAT(I))-1.D0+ALPHA
61      B(I)=DSQRT(DBLE(FLOAT(I))*(DBLE(FLOAT(I))+ALPHA))
	A(N)=2.D0*DBLE(FLOAT(N))-1.D0+ALPHA
70      RETURN
	END
C
	SUBROUTINE IMTQL2(N,D,E,Z,IERR)
C
C       Questa routine implementa il metodo QL e determina gli
C       autovalori e la prima componente dei corrispondenti auto-
C       vettori ortonormali associati alla matrice tridiagonale
C       simmetrica definita dai coefficienti A(I) e B(I) calco-
C       lati con la subroutine CLASS
C
	DOUBLE PRECISION D,E,Z,B,C,F,G,P,R,S,MACHEP,DSQRT,DABS,DSIGN
	DIMENSION D(N),E(N),Z(N)
C
C       MACHEP rappresenta la precisione relativa di macchina ed e' un
C       parametro proprio del calcolatore su cui si lavora. Per esempio,
C       sul Digital PDP 11 abbiamo MACHEP=2.D0**(-56)
C
	DATA MACHEP/1.D-16/
	IERR=0
	IF(N.EQ.1)GO TO 1001
	E(N)=0.D0
	DO 240 L=1,N
	J=0
105     DO 110 M=L,N
	IF(M.EQ.N)GOTO 120
	IF(DABS(E(M)).LE.MACHEP*(DABS(D(M))+DABS(D(M+1))))GOTO 120
110     CONTINUE
120     P=D(L)
	IF(M.EQ.L)GOTO 240
	IF(J.EQ.30)GOTO 1000
	J=J+1
	G=(D(L+1)-P)/(2.D0*E(L))
	R=DSQRT(G*G+1.D0)
	G=D(M)-P+E(L)/(G+DSIGN(R,G))
	S=1.D0
	C=1.D0
	P=0.D0
	MML=M-L
	DO 200 II=1,MML
	I=M-II
	F=S*E(I)
	B=C*E(I)
	IF(DABS(F).LT.DABS(G))GOTO 150
	C=G/F
	R=DSQRT(C*C+1.D0)
	E(I+1)=F*R
	S=1.D0/R
	C=C*S
	GOTO 160
150     S=F/G
	R=DSQRT(S*S+1.D0)
	E(I+1)=G*R
	C=1.D0/R
	S=S*C
160     G=D(I+1)-P
	R=(D(I)-G)*S+2.D0*C*B
	P=S*R
	D(I+1)=G+P
	G=C*R-B
	F=Z(I+1)
	Z(I+1)=S*Z(I)+C*F
200     Z(I)=C*Z(I)-S*F
	D(L)=D(L)-P
	E(L)=G
	E(M)=0.D0
	GOTO 105
240     CONTINUE
	DO 300 II=2,N
	I=II-1
	K=I
	P=D(I)
	DO 260 J=II,N
	IF(D(J).GE.P)GOTO 260
	K=J
	P=D(J)
260     CONTINUE
	IF(K.EQ.I)GOTO 300
	D(K)=D(I)
	D(I)=P
	P=Z(I)
	Z(I)=Z(K)
	Z(K)=P
300     CONTINUE
	GOTO 1001
1000    IERR=L
1001    RETURN
	END
C
C
	DOUBLE PRECISION FUNCTION GAMMA(X)
C
C       Routine di servizio per la valutazione della funzione Gamma
C
	DOUBLE PRECISION X,ALGA,DEXP
	IF(X.LT..5D0)GOTO 10
	GAMMA=DEXP(ALGA(X))
	RETURN
10      GAMMA=DEXP(ALGA(X+1.D0))/X
	RETURN
	END
C
	DOUBLE PRECISION FUNCTION ALGA(X)
	DOUBLE PRECISION CNUM,CDEN,XI,DINT,X,XE,SNUM,SDEN,DLOG,P,DBLE
	DIMENSION CNUM(8),CDEN(8)
	DATA CNUM/4.12084318584777D0,85.68982062831317D0,
     *           243.175243524421D0,-261.7218583856145D0,
     *           -922.2613728801522D0,-517.6383498023218D0,
     *           -77.41064071332953D0,-2.208843997216182D0/
	DATA CDEN/1.D0,45.64677187585908D0,377.8372484823942D0,
     *           951.323597679706D0,846.0755362020782D0,
     *           262.3083470269460D0,24.43519662506312D0,
     *           0.4097792921092615D0/
	XI=DINT(X)
	IF(X-XI.GT..5D0) XI=XI+1.D0
	M=IFIX(SNGL(XI))-1.D0
	XE=X
	IF(M.EQ.-1) XE=X+1.D0
	IF(M.GT.0) XE=X-DBLE(FLOAT(M))
	SNUM=CNUM(1)
	SDEN=CDEN(1)
	DO 10 K=2,8
	SNUM=XE*SNUM+CNUM(K)
10      SDEN=XE*SDEN+CDEN(K)
	ALGA=(XE-1.D0)*SNUM/SDEN
	IF(M.GT.-1) GOTO 20
	ALGA=ALGA-DLOG(X)
	RETURN
20      IF(M.EQ.0) RETURN
	P=XE
	IF(M.EQ.1) GOTO 40
	MM1=M-1
C
C       La parte che segue previene possibili overflow nel calcolo
C       di P. In particolare nella if che segue, del tipo
C       IF(M.GT.KK)..., KK e' un numero intero caratteristico della
C       macchina in uso, tale che P non possa andare in overflow.
C
	IF(M.GT.33) GOTO 50
	DO 30 K=1,MM1
30      P=(XE+DBLE(FLOAT(K)))*P
40      ALGA=ALGA+DLOG(P)
	RETURN
50      ALGA=ALGA+DLOG(XE)
	DO 60 K=1,MM1
60      ALGA=ALGA+DLOG(XE+DBLE(FLOAT(K)))
	RETURN
	END