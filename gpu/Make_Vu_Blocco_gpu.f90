Subroutine Make_Vu_Blocco_gpu(i,j,hk,indice_i,indice_j)  
    !Subroutine che costruisce
    !i blocchetti per la matrice
    !Vu.
! i= elemento i-esimo
! j=   " "    j-esimo

USE variable_2Dgeneral
use cudafor
IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABILI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Input
INTEGER(kind=4),INTENT(IN):: i, j, hk, indice_i, indice_j

!Variabili locali
INTEGER(kind=4):: lf, lg

!!!!!!!!!!!!!!!!!!!!!!!!!! CORPO della SUBROUTINE !!!!!!!!!!!!!!!!!!

!Inizializzazione matrice blocco
DO lg=1,grado_q+1
    DO lf=1,grado_q+1
        call Vuintegra_gpu(i,lf,j,lg,hk,indice_i,indice_j)
    END DO
END DO

RETURN
END SUBROUTINE Make_Vu_Blocco_gpu