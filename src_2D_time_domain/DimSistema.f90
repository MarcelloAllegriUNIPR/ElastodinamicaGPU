SUBROUTINE DimSistema

  USE variable_2Dgeneral

  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABILI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  !Variabili locali
  INTEGER(kind=4):: j, AllocateStatus

  !Variabili locali
  LOGICAL(kind=4),EXTERNAL:: ULeftAss
  !!!!!!!!!!!!!!!!!!!!!!!!!! CORPO della SUBROUTINE !!!!!!!!!!!!!!!!!!

  ! Dimensionamento della matrice Vu -> DimVu
  	DimVu=0

	IF(grado_q.EQ.INT(0,4)) THEN    !funzioni di forma costanti
        DimVu=number_elements
    ELSE                            !funzioni di forma di grado >0
        DO j=1,number_elements
            IF (ULeftAss(j)) THEN
	            DimVu=DimVu+grado_q
	        ELSE
	            DimVu=DimVu+grado_q+1
	        END IF
        END DO
    END IF
  
  !ALLOCATE(Vu(DimVu,DimVu),STAT=AllocateStatus)
  !IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
  !ALLOCATE(matrix(2*DimVu,2*DimVu),STAT=AllocateStatus)     !DA CAMBIARE IN matrix(2*DimVu,2*DimVu)
  !IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
  !ALLOCATE(u_bar(DimVu),STAT=AllocateStatus)
  !IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
  !ALLOCATE(rhs(2*DimVu),STAT=AllocateStatus)          
  !IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
  
RETURN

END SUBROUTINE DimSistema


