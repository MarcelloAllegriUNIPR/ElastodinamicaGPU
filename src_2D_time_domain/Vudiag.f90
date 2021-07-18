      Double precision function Vudiag(li,lj,i,hk,indice_i,indice_j)
!
!     INTEGRALE DELLE FUNZIONI DI FORMA li,lj SULL'ELEMENTO i
!     ELEMENTI SOVRAPPOSTI
!
  USE variable_2Dgeneral
  
  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABILI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  INTEGER(kind=4),INTENT(IN) :: i, li, lj, hk ,indice_i, indice_j
  !REAL(kind=8),INTENT(IN) :: velC, velP
  
  !Variabili locali
  REAL(kind=8),EXTERNAL:: Vudiag_sub, fiu
  REAL(kind=8):: Vudiag_1, Vudiag_2, Vudiag_3, Vudiag_4, estremo
  INTEGER(kind=4) :: ki, Ngauss, AllocateStatus
  REAL(kind=8),DIMENSION(:),ALLOCATABLE:: x, w
  INTEGER(kind=4), DIMENSION(2,2) :: delta_kronecker
  !!!!!!!!!!!!!!!!!!!!!!!!!! CORPO della SUBROUTINE !!!!!!!!!!!!!!!!!!
     
	Vudiag=0.d0

	Vudiag_1=Vudiag_sub(li,lj,i,hk*deltat,REAL(0,8),indice_i,indice_j)
	Vudiag_2=Vudiag_sub(li,lj,i,hk*deltat,deltat,indice_i,indice_j)
	Vudiag_3=Vudiag_sub(li,lj,i,(hk+1)*deltat,REAL(0,8),indice_i,indice_j)
	Vudiag_4=Vudiag_sub(li,lj,i,(hk+1)*deltat,deltat,indice_i,indice_j)

	Vudiag=-(1.d0/(8.d0*datan(1.d0)*rho))*(Vudiag_1-Vudiag_2-Vudiag_3+Vudiag_4) !ATTENZIONE al SEGNO MENO!!!!!!!!!
	
	RETURN
	END 
