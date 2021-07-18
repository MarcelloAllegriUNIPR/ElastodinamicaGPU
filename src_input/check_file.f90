SUBROUTINE check_file(file)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                    
! SUBROUTINE NAME: check_file
!
! GOAL           : verifica l'esistenza di un file  
!                                                                    
! INPUTS :  file : nome del file di cui si vuole verificare l'esistenza
! 
! AUTORE: Luca Desiderio                                             
!                                                                    
! DATA:   17/XI/2017                           
!                                                                    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABILI !!!!!!!!!!!!!!!!!!!!!!!!!!  
  
  !Input
  CHARACTER(len=*),INTENT(IN):: file
  
  !Variabili locali
  LOGICAL:: flok

  !!!!!!!!!!!!!!!!!!!!!!!!!! CORPO della SUBROUTINE !!!!!!!!!!!!!!!!!!  

  INQUIRE(FILE=file,EXIST=flok)
   
  IF(.NOT.flok) THEN
     WRITE(22,'(A)') ''
     WRITE(22,'(A,A)') ' Not able to find: ', file
     WRITE(22,'(A)') ''
     WRITE(22,'(A)') 'ANALYSIS TERMINATED'
     PRINT*,'ERROR, file does not exist'
     STOP
  ENDIF

END SUBROUTINE check_file
