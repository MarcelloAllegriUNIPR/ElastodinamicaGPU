SUBROUTINE nodi_pesi_gauss()

  USE variable_2Dgeneral
  
  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABILI !!!!!!!!!!!!!!!!!!!!!!!!!

  !Variabili locali
  INTEGER(kind=4):: Ngauss, dim_gauss, max_gauss, ind_gauss
  INTEGER(kind=4):: AllocateStatus
  
  REAL(kind=8),DIMENSION(2)::endpts
  REAL(kind=8),DIMENSION(:),ALLOCATABLE:: serv, noditemp, pesitemp 

  !!!!!!!!!!!!!!!!!!!!!!!!!! CORPO della SUBROUTINE !!!!!!!!!!!!!!!!!!

  !Allocazione dell'array gauss: Ngauss=2^0...2^8=256
  dim_gauss=9
  max_gauss=2**(dim_gauss-1) 
  
  ALLOCATE(gauss(dim_gauss),STAT=AllocateStatus)
  IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
  
  ALLOCATE(serv(max_gauss),STAT=AllocateStatus)
  IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
  ALLOCATE(noditemp(max_gauss),STAT=AllocateStatus)
  IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
  ALLOCATE(pesitemp(max_gauss),STAT=AllocateStatus)
  IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
  
  DO ind_gauss=1,dim_gauss
  
        Ngauss=2**(ind_gauss-1)
        ALLOCATE(gauss(ind_gauss)%nodiquad(Ngauss),STAT=AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
        ALLOCATE(gauss(ind_gauss)%pesiquad(Ngauss),STAT=AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
        
        CALL gaussq(1,Ngauss,0.d0,0.d0,0,endpts,serv(1:Ngauss),noditemp(1:Ngauss),pesitemp(1:Ngauss))
        
        gauss(ind_gauss)%nodiquad=noditemp(1:Ngauss)
        gauss(ind_gauss)%pesiquad=pesitemp(1:Ngauss)
        
  END DO

RETURN
  
END SUBROUTINE nodi_pesi_gauss
