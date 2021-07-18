SUBROUTINE time2D_input_node(file_num)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                    
! SUBROUTINE NAME: input_node_2D                             
!
! GOAL   : legge le coordinate dei nodi della mesh da un file con
!          estensione .mesh  
!  
! INPUT  : file_num : indice indicativo del file di input
!
! REMARK : i nodi sono memorizzati nel seguente modo
!          coord_x coord_y coord_z ind_domain
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
  INTEGER(kind=4):: ind_nodes
  INTEGER(kind=4):: AllocateStatus

  !!!!!!!!!!!!!!!!!!!!!!!!!! CORPO della SUBROUTINE !!!!!!!!!!!!!!!!!!

  !Lettura del numero di nodi della mesh
  READ(file_num,*) number_nodes

  !Allocazione dell'array list_nodes
  ALLOCATE(list_nodes(number_nodes),STAT=AllocateStatus)
  IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

  !Lettura delle coordinate dei nodi e dell'indice di dominio
  DO ind_nodes=1,number_nodes
     READ(file_num,*) list_nodes(ind_nodes)%coordinates(1:3),&
          list_nodes(ind_nodes)%ind_domain
  END DO

  RETURN
  
END SUBROUTINE time2D_input_node
