SUBROUTINE time2D_input_edges(file_num)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                    
! SUBROUTINE NAME: input_edges                             
!
! GOAL   : legge le incidenze dei segmenti di bordo della mesh da un 
!          file con estensione .mesh  
!  
! INPUT  : file_num : indice indicativo del file di input
!
! REMARK : i segmenti sono memorizzati nel seguente modo
!          incid_1 incid_2 ind_BD
!     
! AUTORE: Luca Desiderio                                             
!                                                                    
! DATA:   XI/13/2017                           
!                                                                    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

  USE variable_2Dgeneral
  
  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABILI !!!!!!!!!!!!!!!!!!!!!!!!!
  
  !Input
  INTEGER(kind=4),INTENT(IN):: file_num

  !Variabili locali
  INTEGER(kind=4):: ind_elements
  INTEGER(kind=4):: AllocateStatus

  !!!!!!!!!!!!!!!!!!!!!!!!!! CORPO della SUBROUTINE !!!!!!!!!!!!!!!!!!
 
  READ(file_num,*) number_elements

  !Allocazione dell'array list_elements
  ALLOCATE(list_elements(number_elements),STAT=AllocateStatus)
  IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

  !Lettura delle incidenze dei triangoli
  DO ind_elements=1,number_elements
     READ(file_num,*) list_elements(ind_elements)%nodes(1:2),&
          list_elements(ind_elements)%ind_BD
  END DO

  RETURN

END SUBROUTINE time2D_input_edges
