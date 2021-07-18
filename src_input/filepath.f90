SUBROUTINE filepath(path,dir,file,ext)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                    
! SUBROUTINE NAME: filepath
!
! GOAL:  la subroutine crea il path necessario per reperire un file nel
!        file system
!
! INPUT: dir,   nome della directory in cui si trova il file
!        file,  nome del file
!        ext,   estensione del file
!                                                                    
! OUTPUT: path, path del file
! 
! AUTORE: Luca Desiderio                                             
!                                                                    
! DATA:   1/II/2018                         
!                                                                    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
	
  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABILI !!!!!!!!!!!!!!!!!!!!!!!!!!

  !Variabili di Input
  CHARACTER(LEN=*),INTENT(IN):: dir, file, ext

  !Variabili di Output
  CHARACTER(LEN=*),INTENT(OUT):: path

  !!!!!!!!!!!!!!!!!!!!!!!!!! CORPO della SUBROUTINE !!!!!!!!!!!!!!!!!!

  path=TRIM(dir) // '/' // TRIM(file) // TRIM(ext)

  RETURN
 
END SUBROUTINE filepath
