SUBROUTINE time2D_toeplitz_RHS(file_output)

  USE OMP_LIB !abilitazione openmp --> proprietà progetto --> proprietà di config --> fortran --> Language --> Process OpenMP Directives (Generate Parallel Code /QOpenMP)
  USE variable_2Dgeneral
  USE Variables
  use kernels
  USE VuExtraGpu
  use cudafor
  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABILI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  CHARACTER(len=200),INTENT(IN):: file_output
  
  !Variabili locali
  INTEGER(kind=4):: rang, i_time, i, j
  INTEGER(kind=4):: AllocateStatus
  Integer :: thread,tstart, tstop, rate, ind_gauss, istat
  !CHARACTER(len=8):: fmt  !Descrittore di formato
  CHARACTER(len=8):: i_time_st
  
  CHARACTER(len=200)::  file_mat, file_RHS !dir, ext
  CHARACTER(len=200)::  file_path_mat, file_path_RHS
  DOUBLE PRECISION START
  !Calcolo della dimensione della matrice con allocazione di matrix e Vu
  CALL DimSistema
  
  !Formato per la conversione da integer a character della variabile i_time
  !(l'intero viene trasformato in una stringa di lunghezza 5 con 0 a sx)  
  fmt = '(I5.5)'
  ind_gauss = 6
  ext='.txt'


  if (useGpu .eq. 0) then
    dir='./../tests/output/cpu'
  else
    dir='./../tests/output/gpu'
    NGaussDimension = 2**(ind_gauss-1)
    
    dimGrid=dim3(1,1,1)
    dimBlockPreCalculation = dim3(NGaussDimension,1,1)
    dimBlockCalculation = dim3(NGaussDimension,NGaussDimension,1)
    VuExtraGrid = dim3(1,1,1)
    VuExtraBlock= dim3(2*(ind_gauss-1),4,1)
	  
    dimension_d = NGaussDimension

    call setCommonData(gauss(ind_gauss)%pesiquad, gauss(ind_gauss)%pesiquad, &
                       gauss(ind_gauss)%nodiquad, gauss(ind_gauss)%nodiquad, &
                       velC_P, velC_S, grado_q)    
    rho_d = rho
    DimVu_d = DimVu

    do i = 1 , 2*(ind_gauss-1)*4
      istat = cudaStreamCreate(stream(i))
    enddo

    delta_kronecker_VuExtra(1,1)=1  
    delta_kronecker_VuExtra(1,2)=0  
    delta_kronecker_VuExtra(2,1)=0  
    delta_kronecker_VuExtra(2,2)=1

    ALLOCATE(x_VuExtra(2**(ind_gauss-1)),STAT=AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory ***" 
  
    x_VuExtra=gauss(ind_gauss)%nodiquad 

    allocate (xtrasl(NGaussDimension), A1(NGaussDimension), B1(NGaussDimension), alfa_j1(NGaussDimension), beta_j1(NGaussDimension),&
     xinttrasl(NGaussDimension),s(NGaussDimension),ds(NGaussDimension))

     sharedMemDimension = NGaussDimension*NGaussDimension*sizeof(START)
  endif
  

START = omp_get_wtime()
if(useGpu .eq. 0) then
  DO i_time=250,250
    CALL time2D_toeplitz(i_time,file_output);    
  end do
else
   DO i_time=250,250
       CALL time2D_toeplitz_gpu(i_time,file_output);    
   end do
endif

!DO i_time=250,300
!    CALL time2D_RHS(i_time,file_output)    
!END DO

if (useGpu .eq. 0) then
  print *, "cpu ", omp_get_wtime() - START, " seconds."
  print *, "cpu only VuExtra cycles", cputime, " seconds."
else 
  print *, "gpu ", omp_get_wtime() - START, " seconds."
  print *, "gpu only VuExtra cycles", gputime, " seconds."
endif
!pause
!print *, "cpuwins:", cpuwins, "gpuwins:", gpuwins
RETURN
END SUBROUTINE time2D_toeplitz_RHS
