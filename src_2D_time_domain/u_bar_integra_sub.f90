      Double precision function u_bar_integra_sub(i,l,tempo1,indice_termine_noto)
  USE variable_2Dgeneral
  
  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABILI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  INTEGER(kind=4),INTENT(IN):: i, l, indice_termine_noto
  REAL(kind=8),INTENT(IN):: tempo1
  
  !Variabili locali
  INTEGER(kind=4):: ind_gauss, Ngauss, AllocateStatus, ki, DirDat
   
  REAL(kind=8)::len_i, sum
  
  REAL(kind=8),DIMENSION(:),ALLOCATABLE:: x, w
  
  REAL(kind=8),EXTERNAL::fiU, posx, posy, Unota
  
  !!!!!!!!!!!!!!!!!!!!!!!!!! CORPO della SUBROUTINE !!!!!!!!!!!!!!!!!!
	
	ind_gauss=7
    Ngauss=2**(ind_gauss-1) !32 nodi di Gauss
	
    ALLOCATE(x(Ngauss),STAT=AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    ALLOCATE(w(Ngauss),STAT=AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    
    x=gauss(ind_gauss)%nodiquad
    w=gauss(ind_gauss)%pesiquad

	sum=0.d0
	len_i=(list_elements(i)%length)/2.d0
	DirDat=list_elements(i)%ind_BD

    DO ki=1,Ngauss
	   sum=sum+w(ki)*Unota(DirDat,posx(i,x(ki)),posy(i,x(ki)),tempo1,indice_termine_noto)*len_i*fiU(l,len_i*(x(ki)+1.d0),len_i*2.d0,grado_q) !!!,indice_termine_noto) da aggiungere in unota!!!
      !sum=sum+w(ki)*len_i*fiU(l,len_i*(x(ki)+1.d0),len_i*2.d0,grado_q)*Unota(DirDat,posx(i,x(ki)),posy(i,x(ki)),tempo1,indice_termine_noto)
	  !write(*,*) 'Unota, nodo', Unota(DirDat,posx(i,x(ki)),posy(i,x(ki)),tempo1,indice_termine_noto), ki
	END DO
	 ! write(*,*) 'elemento', i
	 ! pause

    u_bar_integra_sub=sum*0.5d0
	
	

    RETURN
    END