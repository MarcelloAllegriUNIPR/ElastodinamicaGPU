SUBROUTINE time2D_solver(file_output)
! Risoluzione del sistema lineare con time marching

  USE variable_2Dgeneral

  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABILI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  CHARACTER(len=200),INTENT(IN):: file_output
  
  !Variabili locali
  INTEGER(kind=4):: i_time, j_time, i, j
  INTEGER(kind=4):: AllocateStatus, INFO
  
  INTEGER(kind=4),DIMENSION(:),ALLOCATABLE:: IPIV
  
  REAL(kind=8),DIMENSION(2*dimVu,2*dimVu):: matrices !!!!!!!!!!!!!!!!!---MIA MODIFICA
  REAL(kind=8),DIMENSION(:,:),ALLOCATABLE:: matrix
  !CHARACTER(len=8):: fmt  !Descrittore di formato
  CHARACTER(len=8):: i_time_st, j_time_st
  
  CHARACTER(len=200)::  file_mat, file_RHS, file_sol !dir, ext
  CHARACTER(len=200)::  file_path_mat, file_path_RHS, file_path_sol
  
  !!!!!!!!!!!!!!!!!!!!!!!!!! CORPO della SUBROUTINE !!!!!!!!!!!!!!!!!!
  
  !Allocazione del vettore soluzione
  ALLOCATE(sol(2*dimVu*Nt),STAT=AllocateStatus)!!!!!!!!!!!!!!!!!---MIA MODIFICA
  IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
  ALLOCATE(matrix(2*DimVu,2*DimVu),STAT=AllocateStatus)
  IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
  !Nome della directory in cui si legge il file contenente il blocco principale
  !dir='./../tests/output'
  !Estensione del file
  !ext='.txt'
 
  !Formato per la conversione da integer a character della variabile i_time
  !(l'intero viene trasformato in una stringa di lunghezza 5 con 0 a sx)
  !fmt = '(I5.5)'
  
  !Inizializzazione della variabile i_time
  i_time=1
  
  !Conversione di i_time da integer a character
  WRITE (i_time_st,fmt) i_time
  !Nome del file
  file_mat=TRIM(file_output) // '_matrice_' // TRIM(i_time_st)
  !Costruzione del path del file contenente il blocco principale ed il relativo RHS
  CALL filepath(file_path_mat,dir,file_mat,ext)
  !Apertura del file contenente il blocco principale
  OPEN(i_time,FILE=TRIM(file_path_mat))
  !Lettura del blocco principale
  DO i=1,2*DimVu
    DO j=1,2*DimVu
        READ(i_time,*) matrix(i,j)
    END DO
  END DO
  CLOSE(i_time)
  
  
  !Allocazione dell'array IPIV
  ALLOCATE(IPIV(2*dimVu),STAT=AllocateStatus)
  IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
 
  !Fattorizzazione LU del blocco principale (http://www.netlib.org/lapack/explore-html/dd/d9a/group__double_g_ecomputational_ga0019443faea08275ca60a734d0593e60.html#ga0019443faea08275ca60a734d0593e60)
  CALL DGETRF(2*DimVu,2*DimVu,matrix,2*DimVu,IPIV,INFO)
  IF (INFO /= 0) STOP "*** Error LU FACTORIZATION ***"
  
  !Nome del file
  file_RHS=TRIM(file_output) // '_rhs_' // TRIM(i_time_st)
  !Costruzione del path per il file di in cui si legge il primo RHS 
  CALL filepath(file_path_RHS,dir,file_RHS,ext)
  !Apertura dei file
  OPEN(i_time,FILE=TRIM(file_path_RHS))
  !Lettura del primo RHS
  DO i=1,2*DimVu
    READ(i_time,*) sol(i)
  END DO
  CLOSE(i_time)
  
  !Risoluzione sistema lineare (http://www.netlib.org/lapack/explore-html/dd/d9a/group__double_g_ecomputational_ga58e332cb1b8ab770270843221a48296d.html#ga58e332cb1b8ab770270843221a48296d)
  CALL DGETRS('N',2*dimVu,1,matrix,2*dimVu,IPIV,sol(1:2*dimVu),2*dimVu,INFO)	
  IF (INFO /= 0) STOP "*** Error SYSTEM ***"
  
  !Nome del file in cui viene scritta la soluzione
  file_sol=TRIM(file_output) // '_soluzione'
  !Costruzione del path del file in cui viene scritta la soluzione
  CALL filepath(file_path_sol,dir,file_sol,ext)
  !Apertura del file in cui viene scritta la soluzione
  OPEN(1,FILE=TRIM(file_path_sol))
  !Scrittura della soluzione
  DO i=1,2*DimVu
    WRITE(1,*) sol(i)
  END DO
  
  DO i_time=2,Nt
  
    !Conversione di i_time da integer a character
    WRITE (i_time_st,fmt) i_time
     
	!***** STEP 1: lettura RHS corrente
    !Nome del file
    file_RHS=TRIM(file_output) // '_rhs_' // TRIM(i_time_st)
    !Costruzione del path per il file di in cui si legge il RHS corrente 
    CALL filepath(file_path_RHS,dir,file_RHS,ext)
    !Apertura dei file
    OPEN(i_time,FILE=TRIM(file_path_RHS))
    !Lettura del RHS corrente
    DO i=1,2*DimVu
        READ(i_time,*) sol((i_time-1)*2*dimVu+i)
    END DO
    CLOSE(i_time)
	 
    !***** STEP 2: lettura blocchi e aggiornamento RHS	
	DO j_time=2,i_time
	    WRITE (j_time_st,fmt) j_time   !!!!!aggiungi come variabile carattere j_time_st
	    file_mat=TRIM(file_output) // '_matrice_' // TRIM(j_time_st)
	    CALL filepath(file_path_mat,dir,file_mat,ext)
	    !Apertura del file contenente il blocco secondario corrente
	    OPEN(j_time,FILE=TRIM(file_path_mat))
		!Lettura del blocco secondario corrente
		DO i=1,2*DimVu
		    DO j=1,2*DimVu  
		        READ(j_time,*) matrices(i,j)  !ridurre le dimensioni di matrices
			END DO
        END DO			
		CLOSE(j_time)
		
		!(http://www.netlib.org/lapack/explore-html/d7/d15/group__double__blas__level2_gadd421a107a488d524859b4a64c1901a9.html#gadd421a107a488d524859b4a64c1901a9)
		!*****aggiornamento RHS
		CALL DGEMV('N',2*dimVu,2*dimVu,REAL(-1,8),matrices(:,:),2*dimVu,sol((i_time-j_time)*2*dimVu+1:(i_time-j_time+1)*2*dimVu),&
        1,REAL(1,8),sol((i_time-1)*2*dimVu+1:i_time*2*dimVu),1)
		
	END DO

       
	!***** STEP 3: risoluzione del sistema corrente
	CALL DGETRS('N',2*dimVu,1,matrix,2*dimVu,IPIV,sol((i_time-1)*2*dimVu+1:i_time*2*dimVu),2*dimVu,INFO)	
    IF (INFO /= 0) STOP "*** Error SYSTEM ***"
   
    !***** STEP 4: scrittura della soluzione corrente
    DO i=1,2*DimVu
        WRITE(1,*) sol((i_time-1)*2*dimVu+i)
    END DO
  
  END DO
  
  CLOSE(1)
  
  RETURN
  
  END SUBROUTINE time2D_solver