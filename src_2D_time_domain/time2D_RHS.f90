SUBROUTINE time2D_RHS(i_time,file_output)

  USE variable_2Dgeneral
  USE OMP_LIB
  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABILI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  INTEGER(kind=4),INTENT(IN):: i_time
  CHARACTER(len=200),INTENT(IN):: file_output
  !Variabili locali
  INTEGER(kind=4):: i, iblocco, iriga, hk,AllocateStatus
  REAL(kind=8),DIMENSION(grado_q+1):: blocco
  LOGICAL(kind=4),EXTERNAL:: ULeftAss
  CHARACTER(len=200)::  file_RHS, file_path_RHS
  CHARACTER(len=8):: i_time_st
  REAL(kind=8),DIMENSION(:),ALLOCATABLE:: u_bar
  REAL(kind=8),DIMENSION(:),ALLOCATABLE:: rhs 
  !!!!!!!!!!!!!!!!!!!!!!!!!! CORPO della SUBROUTINE !!!!!!!!!!!!!!!!!!
  hk=i_time-1
  
  ALLOCATE(u_bar(DimVu),STAT=AllocateStatus)
  IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
  ALLOCATE(rhs(2*DimVu),STAT=AllocateStatus)          
  IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
  
  !--------PARTE 1 TERMINE NOTO----------!
  CALL DSCAL(DimVu,REAL(0,8),u_bar,1)
     
  !Creazione degli elementi del termine noto
  iblocco=grado_q+1
  iriga=1
  DO i=1,number_elements
     !***** STEP 1: Numerazione righe termine noto
	 IF ((grado_q.NE.INT(0,4)).and.ULeftAss(i)) THEN
	    iriga=iriga-1
     END IF
     !***** STEP 2: Creazione dei sottoblocchi
     call Make_u_bar_Blocco(i,blocco,hk,1)
	 !***** STEP 3: Assemblaggio
	 call u_bar_copiaBlocco(iriga,blocco,u_bar)
	 iriga=iriga+iblocco
  END DO

  CALL DSCAL(DimVu,REAL(2,8),u_bar,1) !!!ATTENZIONE: moltiplico per 2.d0
  CALL DCOPY(DimVu,u_bar,1,rhs(1:DimVu),1)
  
  !--------PARTE 2 TERMINE NOTO----------!
  CALL DSCAL(DimVu,REAL(0,8),u_bar,1)
     
  !Creazione degli elementi del termine noto
  iblocco=grado_q+1
  iriga=1
  DO i=1,number_elements
     !***** STEP 1: Numerazione righe termine noto
	 IF ((grado_q.NE.INT(0,4)).and.ULeftAss(i)) THEN
	    iriga=iriga-1
     END IF
     !***** STEP 2: Creazione dei sottoblocchi
     call Make_u_bar_Blocco(i,blocco,hk,2)
	 !***** STEP 3: Assemblaggio
	 call u_bar_copiaBlocco(iriga,blocco,u_bar)
	 iriga=iriga+iblocco
  END DO

  CALL DSCAL(DimVu,REAL(2,8),u_bar,1) !!!ATTENZIONE: moltiplico per 2.d0
  CALL DCOPY(DimVu,u_bar,1,rhs(DimVu+1:2*DimVu),1)
  
  !write(*,*) 'scrittura termine noto', i_time, 'thread', OMP_GET_THREAD_NUM ()
  WRITE (i_time_st,fmt) i_time
  file_RHS=TRIM(file_output) // '_rhs_' // TRIM(i_time_st)
  CALL filepath(file_path_RHS,dir,file_RHS,ext)
  OPEN(i_time+Nt+6,FILE=TRIM(file_path_RHS))    
  DO i=1,2*DimVu
      WRITE(i_time+Nt+6,*) rhs(i),i
  END DO
  CLOSE(i_time+Nt+6)
  
END SUBROUTINE time2D_RHS


