      Subroutine Make_u_bar_Blocco(i,blocco,hk,indice_termine_noto)      
                                            !Subroutine che costruisce
                                            !i blocchetti per la matrice
                                            !Vu.
      ! i= elemento i-esimo
      ! j=   " "    j-esimo
      
      USE variable_2Dgeneral
      
      IMPLICIT NONE

      !Input
      INTEGER(kind=4),INTENT(IN):: i, hk, indice_termine_noto
      
      !Input/Output
      REAL(kind=8),DIMENSION(grado_q+1),INTENT(INOUT):: blocco
      
      !Variabili locali
      INTEGER(kind=4):: lg
      REAL(kind=8),EXTERNAL:: u_bar_integra
    
      !Inizializzazione blocco termine noto
      CALL DSCAL(grado_q+1,0.d0,blocco,1)

      DO lg=1,grado_q+1
         blocco(lg)=-u_bar_integra(i,lg,hk,indice_termine_noto) !!!! ATTENZIONE al SEGNO
		 ! write(*,*) 'elemento', i
		 ! write(*,*) 'grado_q', grado_q
		 ! write(*,*) 'funzione di forma', lg
		 ! write(*,*) 'risultato integrale', blocco(lg)
		 ! pause
      END DO

      RETURN
      END SUBROUTINE Make_u_bar_Blocco
