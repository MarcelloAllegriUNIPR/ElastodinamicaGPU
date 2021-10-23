SUBROUTINE time2D_toeplitz_gpu(i_time,file_output)

	USE variable_2Dgeneral
	use Variables
	use kernels
	use VuExtraGpu
	IMPLICIT NONE
  
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABILI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
	!Input
	INTEGER(kind=4),INTENT(IN):: i_time
	CHARACTER(len=200),INTENT(IN):: file_output
	
	!Variabili locali
	INTEGER(kind=4):: i, j, iblocco, jblocco, iriga, icolonna, hk, indice_i, indice_j,AllocateStatus
	REAL(kind=8),DIMENSION(grado_q+1,grado_q+1):: blocco
	double precision :: app
	
	LOGICAL(kind=4),EXTERNAL:: ULeftAss  
	CHARACTER(len=200)::  file_mat, file_path_mat
	CHARACTER(len=8):: i_time_st
	type (dim3) :: grid, tBlock
    integer :: TILE_DIM = 32, BLOCK_ROWS = 8
	
	!!!!!!!!!!!!!!!!!!!!!!!!!! CORPO della SUBROUTINE !!!!!!!!!!!!!!!!!!
  
	grid = dim3(DimVu/TILE_DIM, DimVu/TILE_DIM, 1)
    tBlock = dim3(TILE_DIM, BLOCK_ROWS, 1)

	ALLOCATE(Vu_device(DimVu,DimVu),STAT=AllocateStatus)
	IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
	ALLOCATE(matrix(2*DimVu,2*DimVu),STAT=AllocateStatus)     !DA CAMBIARE IN matrix(2*DimVu,2*DimVu)
	IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
	
	hk=i_time-1
	delta_x_VuExtra(1) = hk*deltat
	delta_x_VuExtra(2) = (hk*deltat)-deltat
	delta_x_VuExtra(3) = (hk+1)*deltat
	delta_x_VuExtra(4) = ((hk+1)*deltat)-deltat
	delta_x_d = delta_x_VuExtra
	
	DO indice_i=1,2
	   DO indice_j=1,2
		   Vu_device = 0.d0
		   coeff_delta_kronecker_d = delta_kronecker_VuExtra(indice_i,indice_j)/2.d0
		   indice_i_d = indice_i
		   indice_j_d = indice_j
			
		   !Creazione degli elementi della matrice  
		   iblocco=grado_q+1
		   jblocco=grado_q+1
		   iriga=1

		   DO i=1,1                           
				!***** STEP 1: Numerazione righe e colonne matrice
			   IF ((grado_q.NE.INT(0,4)).and.ULeftAss(i)) THEN
				   iriga=iriga-1
			   END IF
  
			   icolonna=1
			   DO j=1,number_elements
					IF ((grado_q.NE.INT(0,4)).and.ULeftAss(j)) THEN
						icolonna=icolonna-1
					END IF
					!***** STEP 2: Creazione dei sottoblocchi: funzioni di forma su elemento i-esimo vs j-esimo                 
					call Make_Vu_Blocco_gpu(i,j,hk,indice_i,indice_j)
					! !***** STEP 3: Assemblaggio
					! call VuCopiaBlocco(iriga,icolonna,blocco,Vu)
					icolonna= icolonna+ jblocco
			   END DO
			   iriga=iriga+iblocco            
		   END DO
	
		   call copySharedMem<<<grid, tBlock,sizeof(app)*(TILE_DIM)*(TILE_DIM)>>>(indice_i, indice_j)
		   !matrix_device --> matrix
	   END DO
	END DO
  
	! WRITE (i_time_st,fmt) i_time
	! file_mat=TRIM(file_output) // '_matrice_' // TRIM(i_time_st)
	! CALL filepath(file_path_mat,dir,file_mat,ext)  
	!   OPEN(i_time+6,FILE=TRIM(file_path_mat))
	!      DO i=1,2*DimVu      
	!        DO j=1,2*DimVu    
	!           !print *, matrix(i,j),i,j 
	!           !pause
	!           !WRITE(i_time+6,*) matrix(i,j),i,j 
	!       END DO
	!      END DO
	!    CLOSE(i_time+6)
  
	  DEALLOCATE(Vu_device)
	  DEALLOCATE(matrix_device)
  
  RETURN
  
  END SUBROUTINE time2D_toeplitz_gpu
  
  
  