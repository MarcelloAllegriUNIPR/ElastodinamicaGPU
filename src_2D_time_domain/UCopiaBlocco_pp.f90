Subroutine UCopiaBlocco_pp(icolonna,blocco,suu)
                                            !Subroutine che copia
                                            !i blocchetti per il post processing.
      USE variable_2Dgeneral
      
      IMPLICIT NONE
      
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABILI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !Input
      INTEGER(kind=4),INTENT(IN):: icolonna
      REAL(kind=8),DIMENSION(grado_q+1),INTENT(IN):: blocco
      
      !Output
      REAL(kind=8),DIMENSION(DimVu),INTENT(OUT):: suu

      
      !Variabili locali
      INTEGER(kind=4):: i, icolumn

  !!!!!!!!!!!!!!!!!!!!!!!!!! CORPO della SUBROUTINE !!!!!!!!!!!!!!!!!!

      DO i=1,grado_q+1
            icolumn=icolonna+i-1
            IF (icolumn.eq.0) icolumn=DimVu       !Questa istruzione serve
                                                  !solo nel caso di frontiera
                                                  !semplice

            suu(icolumn)=suu(icolumn)+blocco(i)
      END DO


RETURN
END Subroutine UCopiaBlocco_pp