      Double precision function Vucont(li,lj,i,j,hk,indice_i,indice_j)
!
!     INTEGRALE DELLE FUNZIONI DI FORMA li,lj SUGLI ELEMENTI i,j
!     ELEMENTI CONTIGUI
!      
  USE variable_2Dgeneral
  
  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABILI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  INTEGER(kind=4),INTENT(IN) :: i, j, li, lj, hk, indice_i, indice_j
  !REAL(kind=8),INTENT(IN) :: velC, velP
  
  !Variabili locali
  REAL(kind=8),EXTERNAL:: Vucont_sub
  REAL(kind=8):: Vucont_1, Vucont_2, Vucont_3, Vucont_4

  !!!!!!!!!!!!!!!!!!!!!!!!!! CORPO della SUBROUTINE !!!!!!!!!!!!!!!!!!
    Vucont=0.d0
    
	! If(j.eq.15.and.i.eq.16)then
	! write(*,*) 'elementi', j, i
	! write(*,*) 'indici', indice_i, indice_j
	! write(*,*) 'tempo', hk*deltat
	Vucont_1=Vucont_sub(li,lj,i,j,hk*deltat,REAL(0,8),indice_i,indice_j)
	Vucont_2=Vucont_sub(li,lj,i,j,hk*deltat,deltat,indice_i,indice_j)
	Vucont_3=Vucont_sub(li,lj,i,j,(hk+1)*deltat,REAL(0,8),indice_i,indice_j)
	Vucont_4=Vucont_sub(li,lj,i,j,(hk+1)*deltat,deltat,indice_i,indice_j)
    !Vucont=Vucont_1
	! write(*,*) 'risultato_int', Vucont_1
	! pause
	! endif
	
	
	Vucont=-(1.d0/(8.d0*datan(1.d0)*rho))*(Vucont_1-Vucont_2-Vucont_3+Vucont_4)!ATTENZIONE al SEGNO MENO!!!!!!!!!

    
	RETURN
	END