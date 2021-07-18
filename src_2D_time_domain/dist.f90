SUBROUTINE dist(xest,yest,ind_el)
! Valuta l'elemento ind_el del contorno di minima distanza val_dist dal punto (xest,yest)

  USE variable_2Dgeneral

  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABILI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  REAL(kind=8),INTENT(IN):: xest, yest
  
  !Output
  INTEGER(kind=4),INTENT(OUT):: ind_el
  
  !Variabili locali
  INTEGER(kind=4):: i_elements, ind_nodo1, ind_nodo2
  
  REAL(kind=8):: x_nodo1, y_nodo1, x_nodo2, y_nodo2

  !!!!!!!!!!!!!!!!!!!!!!!!!! CORPO della SUBROUTINE !!!!!!!!!!!!!!!!!!

  ind_el=0
  DO i_elements=1,number_elements
     !Indici degli estremi dell'elemento i_elements 
     ind_nodo1=list_elements(i_elements)%nodes(1)
     ind_nodo2=list_elements(i_elements)%nodes(2)
     !Coordinate degli estremi dell'elemento i_elements 
     x_nodo1=list_nodes(ind_nodo1)%coordinates(1)
     y_nodo1=list_nodes(ind_nodo1)%coordinates(2)
     x_nodo2=list_nodes(ind_nodo2)%coordinates(1)
     y_nodo2=list_nodes(ind_nodo2)%coordinates(2)
     IF((y_nodo1.EQ.y_nodo2).AND.(y_nodo1.EQ.yest).AND.&
     (xest.GE.x_nodo1).AND.(xest.LE.x_nodo2)) THEN
        ind_el=i_elements
        RETURN
     END IF   
     IF((x_nodo1.EQ.x_nodo2).AND.(x_nodo1.EQ.xest).AND.&
     (yest.GE.y_nodo1).AND.(yest.LE.y_nodo2)) THEN
        ind_el=i_elements
        RETURN
     END IF 
     IF((DABS((xest-x_nodo1)/(x_nodo2-x_nodo1)-(yest-y_nodo1)/(y_nodo2-y_nodo1)).LE.1.d-12).AND.&
     (xest.GE.x_nodo1).AND.(xest.LE.x_nodo2)) THEN
        ind_el=i_elements
        RETURN
     END IF
  END DO

  RETURN
  
  END SUBROUTINE dist