      Subroutine Make_Vu_Blocco(i,j,blocco,hk,indice_i,indice_j)  
                                            !Subroutine che costruisce
                                            !i blocchetti per la matrice
                                            !Vu.
      ! i= elemento i-esimo
      ! j=   " "    j-esimo
      
      USE variable_2Dgeneral
      
      IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABILI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !Input
      INTEGER(kind=4),INTENT(IN):: i, j, hk, indice_i, indice_j
      
      !Input/Output
      REAL(kind=8),DIMENSION(grado_q+1,grado_q+1),INTENT(OUT):: blocco
      
      !Variabili locali
      INTEGER(kind=4):: lf, lg
      REAL(kind=8),EXTERNAL:: Vuintegra
    
  !!!!!!!!!!!!!!!!!!!!!!!!!! CORPO della SUBROUTINE !!!!!!!!!!!!!!!!!!

      !Inizializzazione matrice blocco
      DO lg=1,grado_q+1
        CALL DSCAL(grado_q+1,0.d0,blocco(:,lg),1)
      END DO

      DO lg=1,grado_q+1
         DO lf=1,grado_q+1
         !   blocco(lf,lg)=Vuintegra(i,lf,j,lg,hk,velC_S,velC_P,indice_i,indice_j)
         blocco(lf,lg)=Vuintegra(i,lf,j,lg,hk,indice_i,indice_j)
         END DO
      END DO

      RETURN
      END SUBROUTINE Make_Vu_Blocco
