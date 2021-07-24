SUBROUTINE time2D_toeplitz_RHS(file_output)

  USE OMP_LIB !abilitazione openmp --> proprietà progetto --> proprietà di config --> fortran --> Language --> Process OpenMP Directives (Generate Parallel Code /QOpenMP)
  USE variable_2Dgeneral
  USE VuExtraGpu

  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABILI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  CHARACTER(len=200),INTENT(IN):: file_output
  
  !Variabili locali
  INTEGER(kind=4):: rang, i_time, i, j
  INTEGER(kind=4):: AllocateStatus
  Integer :: thread,tstart, tstop, rate, ind_gauss
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
  if (useGpu .eq. 0) then
  dir='./../tests/output/cpu'
  else
    dir='./../tests/output/gpu'
    call setCommonData(gauss(ind_gauss)%pesiquad, gauss(ind_gauss)%pesiquad, &
                       gauss(ind_gauss)%nodiquad, gauss(ind_gauss)%nodiquad, &
                       velC_P, velC_S, grado_q)
  endif
  ext='.txt'

START = omp_get_wtime()
DO i_time=250,250
    CALL time2D_toeplitz(i_time,file_output);    
end do
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
