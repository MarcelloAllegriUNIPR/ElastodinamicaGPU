      Subroutine u_bar_copiaBlocco(iriga,blocco,u_bar) 
                                            !Subroutine che copia
                                            !i blocchetti per il termine noto.
      USE variable_2Dgeneral
      
      IMPLICIT NONE
      
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABILI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !Input
      INTEGER(kind=4),INTENT(IN):: iriga
      REAL(kind=8),DIMENSION(grado_q+1),INTENT(IN):: blocco
      REAL(kind=8),DIMENSION(DimVu),INTENT(INOUT):: u_bar
      
      !Variabili locali
      INTEGER(kind=4):: i, irow

  !!!!!!!!!!!!!!!!!!!!!!!!!! CORPO della SUBROUTINE !!!!!!!!!!!!!!!!!!

      DO i=1,grado_q+1
            irow=iriga+i-1
            IF (irow.eq.0) irow=DimVu             !Questa istruzione serve
                                                  !solo nel caso di frontiera
                                                  !semplice

            u_bar(irow)=u_bar(irow)+blocco(i)
      END DO

      RETURN
      END SUBROUTINE u_bar_copiaBlocco
