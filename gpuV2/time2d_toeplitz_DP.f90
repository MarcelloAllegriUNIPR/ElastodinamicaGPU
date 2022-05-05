SUBROUTINE time2D_toeplitz_DP(i_time,file_output)

    use cudafor
    USE variable_2Dgeneral
    use kernel
    !USE OMP_LIB
    IMPLICIT NONE

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABILI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Input
    INTEGER(kind=4),INTENT(IN):: i_time
    CHARACTER(len=200),INTENT(IN):: file_output
    !Variabili locali
    INTEGER(kind=4):: i, j, iblocco, jblocco, iriga, icolonna, hk, indice_i, indice_j,AllocateStatus,istat,sizeofShMem
    REAL(kind=8),DIMENSION(grado_q+1,grado_q+1):: blocco
    double precision :: var
    REAL(kind=8),DIMENSION(:,:),ALLOCATABLE:: Vu
    REAL(kind=8),DIMENSION(:,:),ALLOCATABLE:: matrix,matrix_cpu
    LOGICAL(kind=4),EXTERNAL:: ULeftAss  
    CHARACTER(len=200)::  file_mat, file_path_mat
    CHARACTER(len=8):: i_time_st
    type (dim3) :: grid, tBlock
    !double precision, device :: hk_d
    !!!!!!!!!!!!!!!!!!!!!!!!!! CORPO della SUBROUTINE !!!!!!!!!!!!!!!!!!

    ALLOCATE(Vu(DimVu,DimVu),STAT=AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    ALLOCATE(matrix(2*DimVu,2*DimVu),STAT=AllocateStatus)     !DA CAMBIARE IN matrix(2*DimVu,2*DimVu)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    ALLOCATE(matrix_cpu(2*DimVu,2*DimVu),STAT=AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
  
    hk=i_time-1
    hk_d = i_time-1
    grid = dim3(2*DimVu,2*DimVu,1)
    !grid = dim3(1,DimVu,1)
    tBlock = dim3(N_gauss,N_gauss/4,1)
    
    sizeofShMem = sizeof(var)*(N_gauss+1)*N_gauss/4
    !call cudaProfilerStart();
    call Make_Vu_Blocco_DP<<<grid,tblock,sizeofShMem>>>(hk_d) 

    ! DO indice_i=1,2
    !   DO indice_j=1,2
    !       DO j=1,DimVu
    !           CALL DSCAL(DimVu,REAL(0,8),Vu(:,j),1)
    !       END DO
 
    !       !Creazione degli elementi della matrice  
    !       iblocco=grado_q+1
    !       jblocco=grado_q+1
    !       iriga=1
    !       DO i=1,number_elements
    !           !call system_clock(tstart, rate);
    !       !***** STEP 1: Numerazione righe e colonne matrice
    !           IF ((grado_q.NE.INT(0,4)).and.ULeftAss(i)) THEN
    !               iriga=iriga-1
    !           END IF

    !           icolonna=1
    !           DO j=1,number_elements
    !               IF ((grado_q.NE.INT(0,4)).and.ULeftAss(j)) THEN
    !                   icolonna=icolonna-1
    !               END IF
    !               !***** STEP 2: Creazione dei sottoblocchi: funzioni di forma su elemento i-esimo vs j-esimo                 
    !               call Make_Vu_Blocco(i,j,blocco,hk,indice_i,indice_j)
                        
    !               !***** STEP 3: Assemblaggio                                    
    !               call VuCopiaBlocco(iriga,icolonna,blocco,Vu)
                
    !               icolonna= icolonna+ jblocco
    !           END DO
    !           iriga=iriga+iblocco            
    !       END DO
 
    !       DO j=1,DimVu
    !           CALL DCOPY(DimVu,Vu(:,j),1,matrix((indice_i-1)*DimVu+1:indice_i*DimVu,(indice_j-1)*DimVu+j),1)    
    !       END DO
     
    !   END DO
    ! END DO

    istat= cudaGetLastError()
    if(istat /= cudaSuccess) write(*,*) 'Sync kernel error: ' , cudaGetErrorString(istat)
    istat= cudaDeviceSynchronize()    
    if(istat /= cudaSuccess) write(*,*) 'Async kernel error : ', cudaGetErrorString(istat)    

    matrix_cpu = matrix_d    
    !istat = cudaDeviceReset()
    !call cudaProfilerStop();

    WRITE (i_time_st,fmt) i_time
  
    ! file_mat=TRIM(file_output) // '_matrice_' // TRIM(i_time_st)
    ! CALL filepath(file_path_mat,dir,file_mat,ext)  
    ! CALL filepath(file_path_mat,dir,file_mat,ext)  
    ! CALL filepath(file_path_mat,dir,file_mat,ext)  
    !  OPEN(i_time+6,FILE=TRIM(file_path_mat))
    !     DO i=1,2*DimVu      
    !       DO j=1,2*DimVu    
    !          !print *, matrix(i,j),i,j 
    !          !pause
    !         if(matrix(i,j) .ne. 0.0d0) then
    !           WRITE(i_time+6,*) matrix(i,j),i,j 
    !         else 
    !           WRITE(i_time+6,*) matrix_cpu(i,j),i,j 
    !         end if
    !      END DO
    !     END DO
    !   CLOSE(i_time+6)
    
    DEALLOCATE(Vu)
    DEALLOCATE(matrix)
    DEALLOCATE(matrix_cpu)
 RETURN
  
END SUBROUTINE time2D_toeplitz_DP
  
  
  