      Double precision function Vuextra(li,lj,i,j,hk,indice_i,indice_j)
!
!     INTEGRALE DELLE FUNZIONI DI FORMA li,lj SUGLI ELEMENTI i,j
!     ELEMENTI SEPARATI
!
  USE variable_2Dgeneral
  
  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABILI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  INTEGER(kind=4),INTENT(IN):: i, j, li, lj, hk, indice_i, indice_j
  !REAL(kind=8),INTENT(IN):: velC, velP
  
  !Variabili locali
  REAL(kind=8),EXTERNAL:: Vuextra_sub
  REAL(kind=8):: Vuextra_1, Vuextra_2, Vuextra_3, Vuextra_4, somme1, somme2, somme

  !!!!!!!!!!!!!!!!!!!!!!!!!! CORPO della SUBROUTINE !!!!!!!!!!!!!!!!!!

	Vuextra=0.d0
    
	Vuextra_1=Vuextra_sub(li,lj,i,j,hk*deltat,REAL(0,8),indice_i,indice_j)
	Vuextra_2=Vuextra_sub(li,lj,i,j,hk*deltat,deltat,indice_i,indice_j)
	Vuextra_3=Vuextra_sub(li,lj,i,j,(hk+1)*deltat,REAL(0,8),indice_i,indice_j)
	Vuextra_4=Vuextra_sub(li,lj,i,j,(hk+1)*deltat,deltat,indice_i,indice_j)  
	
	Vuextra=-(1.d0/(8.d0*datan(1.d0)*rho))*(Vuextra_1-Vuextra_2-Vuextra_3+Vuextra_4) !ATTENZIONE al SEGNO MENO!!!!!!!!!

	RETURN
	END