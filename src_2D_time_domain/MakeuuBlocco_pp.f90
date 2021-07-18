Subroutine MakeuuBlocco_pp(j,blocco,x,y,t,velC,velP,tau,indice_termine_noto,indice_j)    
                                            !Subroutine che costruisce
                                            !i blocchetti per la matrice
                                            !Vu.
      ! i= elemento i-esimo
      ! j=   " "    j-esimo
      
      USE variable_2Dgeneral
      
      IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABILI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !Input
      INTEGER(kind=4),INTENT(IN):: j, indice_termine_noto, indice_j,tau
      REAL(kind=8),INTENT(IN)::x, y, t, velC, velP

     
      !Output
      REAL(kind=8),DIMENSION(grado_q+1),INTENT(OUT):: blocco
      
      !Variabili locali
      INTEGER(kind=4):: lg
      REAL(kind=8),EXTERNAL:: uuextra_pp
    
  !!!!!!!!!!!!!!!!!!!!!!!!!! CORPO della SUBROUTINE !!!!!!!!!!!!!!!!!!

  !Inizializzazione matrice blocco
  CALL DSCAL(grado_q+1,REAL(0,8),blocco,1)
	
  DO lg=1,grado_q+1
    blocco(lg)=uuextra_pp(j,lg,x,y,t,velC,velP,tau,indice_termine_noto,indice_j)
  END DO

RETURN
END SUBROUTINE MakeuuBlocco_pp


