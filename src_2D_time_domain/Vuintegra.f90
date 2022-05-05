      Double precision function Vuintegra(i,li,j,lj,hk,indice_i,indice_j)
!
!     FUNZIONE CHE SI OCCUPA DI INTEGRARE LE FUNZIONI DI FORMA
!     SUI DIVERSI ELEMENTI AL CONTORNO DISTINGUENDO SE SONO
!     ELEMENTI UGUALI, CONTIGUI O SEPARATI
!

  USE variable_2Dgeneral
  
  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABILI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  INTEGER(kind=4),INTENT(IN):: i, li, j, lj, hk, indice_i, indice_j
  !REAL(kind=8),INTENT(IN):: velC, velP

  !Variabili locali
  REAL(kind=8),EXTERNAL:: Vudiag, Vuextra, Vucont

  !!!!!!!!!!!!!!!!!!!!!!!!!! CORPO della SUBROUTINE !!!!!!!!!!!!!!!!!!
      
      IF (i.eq.j) THEN
        Vuintegra= Vudiag(li,lj,j,hk,indice_i,indice_j) !ELASTODINAMICA NORMALE
      ELSE IF(list_elements(i)%nodes(1).EQ.list_elements(j)%nodes(2)) THEN
        Vuintegra= Vucont(li,lj,i,j,hk,indice_i,indice_j) !ELASTODINAMICA NORMALE
      ELSE IF(list_elements(j)%nodes(1).EQ.list_elements(i)%nodes(2)) THEN      
        Vuintegra= Vucont(li,lj,i,j,hk,indice_i,indice_j) !ELASTODINAMICA NORMALE
      ELSE IF((i.eq.number_elements.and.j.eq.1).or.(j.eq.number_elements.and.i.eq.1))then
        Vuintegra= Vucont(li,lj,i,j,hk,indice_i,indice_j) !!!!SOLO PER ARCO CHIUSO
	    ELSE
	      Vuintegra= 0.0d0!Vuextra(li,lj,i,j,hk,indice_i,indice_j) !ELASTODINAMICA NORMALE
	    END IF
	  
 RETURN
 END