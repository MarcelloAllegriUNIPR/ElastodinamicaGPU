      Double precision function posx(k,x)
!
!     FUNZIONE CHE SI OCCUPA DI CALCOLARE LA COORDINATA X  
!     SULL'ELEMENTO DI CONTORNO k A PARTIRE DALLA COORDINATA x
!     SULL'ELEMENTO PARENTE
!

      USE variable_2Dgeneral
  
      IMPLICIT NONE

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABILI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !Input
      INTEGER(kind=4),INTENT(IN):: k
      REAL(kind=8),INTENT(IN):: x
  
      !Variabili locali
      INTEGER(kind=4),DIMENSION(2):: nodi_k
    
      REAL(kind=8),DIMENSION(2):: coor_k1, coor_k2
      
      !!!!!!!!!!!!!!!!!!!!!!!!!! CORPO della SUBROUTINE !!!!!!!!!!!!!!!!!!
      
      !Estremi sottointervallo k 
      nodi_k=list_elements(k)%nodes

      !Coordinate degli estremi 
      coor_k1=list_nodes(nodi_k(1))%coordinates(1:2)
      coor_k2=list_nodes(nodi_k(2))%coordinates(1:2)

      posx = (coor_k1(1)+coor_k2(1))*0.5d0+x*(coor_k2(1)-coor_k1(1))*0.5d0

      RETURN
      END