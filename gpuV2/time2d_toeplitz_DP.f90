SUBROUTINE time2D_toeplitz_DP(i_time,file_output)

    USE variable_2Dgeneral
    use kernel_MakeVuBlocco
    !USE OMP_LIB
    IMPLICIT NONE

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABILI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Input
    INTEGER(kind=4),INTENT(IN):: i_time
    CHARACTER(len=200),INTENT(IN):: file_output
    !Variabili locali
    INTEGER(kind=4):: i, j, iblocco, jblocco, iriga, icolonna, hk, indice_i, indice_j,AllocateStatus
    !REAL(kind=8),DIMENSION(grado_q+1,grado_q+1):: blocco
    double precision :: var
    REAL(kind=8),DIMENSION(:,:),ALLOCATABLE:: Vu
    REAL(kind=8),DIMENSION(:,:),ALLOCATABLE:: matrix
    LOGICAL(kind=4),EXTERNAL:: ULeftAss  
    CHARACTER(len=200)::  file_mat, file_path_mat
    CHARACTER(len=8):: i_time_st
    type (dim3) :: grid, tBlock
    !double precision, device :: hk_d
    !!!!!!!!!!!!!!!!!!!!!!!!!! CORPO della SUBROUTINE !!!!!!!!!!!!!!!!!!

    ! ALLOCATE(Vu(DimVu,DimVu),STAT=AllocateStatus)
    ! IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    ALLOCATE(matrix(2*DimVu,2*DimVu),STAT=AllocateStatus)     !DA CAMBIARE IN matrix(2*DimVu,2*DimVu)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
  
    hk=i_time-1
    hk_d = i_time-1
    !grid = dim3(2*DimVu,2*DimVu,1)
    grid = dim3(1,10,1)
    tBlock = dim3(N_gauss,8,1)
    
    call Make_Vu_Blocco_DP<<<grid,tblock,(sizeof(var)*N_gauss*8)>>>(hk_d) !+sizeof(var)*5+sizeof(hk)*4

    matrix = matrix_d
    print *, matrix(1,1:10)
    !DEALLOCATE(Vu_d)
    !DEALLOCATE(matrix_d)

RETURN
  
END SUBROUTINE time2D_toeplitz_DP
  
  
  