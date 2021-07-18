SUBROUTINE time2D_elabora_input

  USE variable_2Dgeneral
  
  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABILI !!!!!!!!!!!!!!!!!!!!!!!!!!

  !Variabili locali della subroutine
  INTEGER(kind=4):: ind, ind_elements

  INTEGER(kind=4),DIMENSION(2):: nodes

  REAL(kind=8),DIMENSION(3,3):: TF

  !!!!!!!!!!!!!!!!!!!!!!!!!! CORPO della SUBROUTINE !!!!!!!!!!!!!!!!!!

  DO ind_elements=1,number_elements
     nodes=list_elements(ind_elements)%nodes
     DO ind=1,2
        TF(ind,:)=list_nodes(nodes(ind))%coordinates
     END DO
     !Calcolo delle componenti del vettore normale ad ogni elemento
     list_elements(ind_elements)%normal(1)=-(TF(2,2)-TF(1,2))
     list_elements(ind_elements)%normal(2)=TF(2,1)-TF(1,1)
     list_elements(ind_elements)%normal(3)=REAL(0,8)
     !Calcolo della lunghezza di ogni elemento
     list_elements(ind_elements)%length=&
          DSQRT(list_elements(ind_elements)%normal(2)**2+&
          list_elements(ind_elements)%normal(1)**2)
     !Calcolo del versore normale ad ogni elemento
     list_elements(ind_elements)%normal=&
          list_elements(ind_elements)%normal/&
          list_elements(ind_elements)%length
     !Determinazione del tipo di dato noto sull'elemento 
     IF((list_elements(ind_elements)%ind_BD.GE.1).AND.&
          (list_elements(ind_elements)%ind_BD.LE.99)) THEN
        list_elements(ind_elements)%type='u' !!!!!!dato di DIRICHLET
     ELSE IF((list_elements(ind_elements)%ind_BD.GE.101).AND.&
          (list_elements(ind_elements)%ind_BD.LE.199)) THEN
        list_elements(ind_elements)%type='q' !!!!!!dato di NEUMANN
     ELSE
        WRITE(*,*) 'ERROR: invalid boundary data'
        STOP
     END IF
     
  END DO

  RETURN

END SUBROUTINE time2D_elabora_input
