      Double precision function u_bar_integra(i,li,hk,indice_termine_noto)
!
!     FUNZIONE CHE SI OCCUPA DI INTEGRARE LE FUNZIONI DI FORMA
!     SUI DIVERSI ELEMENTI AL CONTORNO DISTINGUENDO SE SONO
!     ELEMENTI UGUALI, CONTIGUI O SEPARATI
!

  USE variable_2Dgeneral
  
  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABILI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  INTEGER(kind=4),INTENT(IN):: i, li, hk, indice_termine_noto
  
  !Variabili locali
  REAL(kind=8),EXTERNAL:: u_bar_integra_sub
  REAL(kind=8)::u_bar_integra_1, u_bar_integra_2

  !!!!!!!!!!!!!!!!!!!!!!!!!! CORPO della SUBROUTINE !!!!!!!!!!!!!!!!!!
      
u_bar_integra_1=u_bar_integra_sub(i,li,hk*deltat,indice_termine_noto)
u_bar_integra_2=u_bar_integra_sub(i,li,(hk+1)*deltat,indice_termine_noto)
      
	  
u_bar_integra=u_bar_integra_1-u_bar_integra_2
	  
 RETURN
 END