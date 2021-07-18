  Double precision function uuextra_pp(j,lj,x,y,t,velC,velP,tau,indice_termine_noto,indice_j)
!
!     INTEGRALE CON FUNZIONE DI FORMA lj SULL'ELEMENTO j
!
  USE variable_2Dgeneral
  
  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABILI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  INTEGER(kind=4),INTENT(IN):: j, lj, tau,indice_termine_noto, indice_j
  REAL(kind=8),INTENT(IN)::velC, velP, x, y, t
  
  !Variabili locali
  REAL(kind=8),EXTERNAL:: uuextra_pp_sub
  REAL(kind=8):: uuextra_1, uuextra_2

  !!!!!!!!!!!!!!!!!!!!!!!!!! CORPO della SUBROUTINE !!!!!!!!!!!!!!!!!!

	uuextra_pp=0.d0

	uuextra_1=uuextra_pp_sub(j,lj,x,y,t,velC,velP,tau*deltat,indice_termine_noto,indice_j)
	uuextra_2=uuextra_pp_sub(j,lj,x,y,t,velC,velP,(tau+1.d0)*deltat,indice_termine_noto,indice_j)
	uuextra_pp=uuextra_1-uuextra_2
    
  RETURN
  END
