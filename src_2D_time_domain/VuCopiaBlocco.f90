      Subroutine VuCopiaBlocco(iriga,icolonna,blocco,Vu) 
                                            !Subroutine che copia
                                            !i blocchetti per la matrice
                                            !Vu.
      USE variable_2Dgeneral
      
      IMPLICIT NONE
      
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABILI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !Input
      INTEGER(kind=4),INTENT(IN):: iriga, icolonna
      REAL(kind=8),DIMENSION(grado_q+1,grado_q+1),INTENT(IN):: blocco
      REAL(kind=8),DIMENSION(DimVu,DimVu),INTENT(INOUT):: Vu
      !Variabili locali
      INTEGER(kind=4):: i, j, irow, icolumn

  !!!!!!!!!!!!!!!!!!!!!!!!!! CORPO della SUBROUTINE !!!!!!!!!!!!!!!!!!
      
      DO i=1,grado_q+1
        DO j=1,grado_q+1
            irow=iriga+i-1
            IF (irow.eq.0) irow=DimVu            !Questa istruzione serve
                                                  !solo nel caso di frontiera
                                                  !semplice
            icolumn=icolonna+j-1
            IF (icolumn.eq.0) icolumn=DimVu          !Questa istruzione serve
                                                  !solo nel caso di frontiera
                                                  !semplice

            Vu(irow,icolumn)=Vu(irow,icolumn)+blocco(i,j) 
        END DO   
      END DO

      RETURN
      END SUBROUTINE VuCopiaBlocco
