subroutine Vuintegra_gpu(i,li,j,lj,hk,indice_i,indice_j)
    !
    !     FUNZIONE CHE SI OCCUPA DI INTEGRARE LE FUNZIONI DI FORMA
    !     SUI DIVERSI ELEMENTI AL CONTORNO DISTINGUENDO SE SONO
    !     ELEMENTI UGUALI, CONTIGUI O SEPARATI
    !
    
      USE variable_2Dgeneral
      use Variables
      IMPLICIT NONE
    
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABILI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REAL(kind=8),EXTERNAL:: Vudiag, Vuextra, Vucont
      !Input
      INTEGER(kind=4),INTENT(IN):: i, li, j, lj, hk, indice_i, indice_j
      double precision :: value
      !!!!!!!!!!!!!!!!!!!!!!!!!! CORPO della SUBROUTINE !!!!!!!!!!!!!!!!!!
          
      value = 0.d0
      
      IF (i.eq.j) THEN
        value = Vudiag(li,lj,j,hk,indice_i,indice_j) !ELASTODINAMICA NORMALE             
      ELSE IF(list_elements(i)%nodes(1).EQ.list_elements(j)%nodes(2)) THEN
        value = Vucont(li,lj,i,j,hk,indice_i,indice_j) !ELASTODINAMICA NORMALE            
      ELSE IF(list_elements(j)%nodes(1).EQ.list_elements(i)%nodes(2)) THEN      
        value = Vucont(li,lj,i,j,hk,indice_i,indice_j) !ELASTODINAMICA NORMALE
      ELSE IF((i.eq.number_elements.and.j.eq.1).or.(j.eq.number_elements.and.i.eq.1)) then
        value = Vucont(li,lj,i,j,hk,indice_i,indice_j) !!!!SOLO PER ARCO CHIUSO
      ELSE
        call Vuextra_gpu(li,lj,i,j,hk,indice_i,indice_j) !ELASTODINAMICA NORMALE
      END IF
      
      if (value .ne. 0.d0) Vu_device(i,j) = value
     RETURN
     END