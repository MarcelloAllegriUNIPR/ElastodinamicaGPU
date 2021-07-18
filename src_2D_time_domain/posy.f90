      Double precision function posy(k,y)
!
!     FUNZIONE CHE SI OCCUPA DI CALCOLARE LA COORDINATA Y  
!     SULL' ELEMENTO AL CONTORNO k A PARTIRE DALLA COORDINATA y
!     SULL'ELEMENTO PARENTE
!

      USE variable_2Dgeneral
  
      IMPLICIT NONE

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABILI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !Input
      INTEGER(kind=4),INTENT(IN):: k
      REAL(kind=8),INTENT(IN):: y
  
      !Variabili locali
      INTEGER(kind=4),DIMENSION(2):: nodi_k
    
      REAL(kind=8),DIMENSION(2):: coor_k1, coor_k2
      
      !!!!!!!!!!!!!!!!!!!!!!!!!! CORPO della SUBROUTINE !!!!!!!!!!!!!!!!!!
      
      !Estremi sottointervallo k 
      nodi_k=list_elements(k)%nodes

      !Coordinate degli estremi 
      coor_k1=list_nodes(nodi_k(1))%coordinates(1:2)
      coor_k2=list_nodes(nodi_k(2))%coordinates(1:2)

      posy = (coor_k1(2)+coor_k2(2))*0.5d0+y*(coor_k2(2)-coor_k1(2))*0.5d0

      RETURN
      END