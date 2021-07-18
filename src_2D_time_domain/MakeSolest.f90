SUBROUTINE MakeSolest(suu,x,y,t,velC,velP,tau,indice_termine_noto,indice_j)    ! Codice per calcolare le matrici Ruu
					
  USE variable_2Dgeneral

  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABILI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  INTEGER(kind=4),INTENT(IN)::tau, indice_termine_noto, indice_j
  REAL(kind=8),INTENT(IN)::x, y, t, velC, velP
  
  !Output
  REAL(kind=8),DIMENSION(DimVu),INTENT(OUT)::suu

  !Variabili locali
  INTEGER(kind=4):: icolonna, jblocco, j
  REAL(kind=8),DIMENSION(grado_q+1)::blocco
  logical(kind=4),EXTERNAL::ULeftAss
  
  !!!!!!!!!!!!!!!!!!!!!!!!!! CORPO della SUBROUTINE !!!!!!!!!!!!!!!!!!

CALL DSCAL(DimVu,REAL(0,8),suu,1)

icolonna=1
jblocco=grado_q+1

DO j=1,number_elements
    !***** STEP 1: Numerazione contributi
    IF ((grado_q.NE.INT(0,4)).and.ULeftAss(j)) THEN
	   icolonna=icolonna-1
    END IF
    !***** STEP 2: Creazione dei sottoblocchi: contributo funzione di forma su elemento j-esimo     
    call MakeuuBlocco_pp(j,blocco,x,y,t,velC,velP,tau,indice_termine_noto,indice_j)
	!***** STEP 3: Assemblaggio
	IF (grado_q.EQ.INT(0,4)) THEN
	    suu(j)=blocco(1)
	else
	    call UCopiaBlocco_pp(icolonna,blocco,suu)
        icolonna= icolonna+ jblocco
	endif
END DO
	
RETURN
END SUBROUTINE MakeSolest

