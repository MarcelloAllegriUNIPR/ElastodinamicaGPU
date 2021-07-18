!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!GOAL:                                                                       !
!                                                                            !
!                                                                            !
!									                                         !
!Luca DESIDERIO & Chiara GUARDASONI							                 !
!luca.desiderio@unipr.it & chiara.guardasoni@unipr.it			             !
!CREATED: 06/03/2018	                   				                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM MAIN

  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABILI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
  !Local Variables
  !INTEGER(kind=4):: nb_threads
  
  CHARACTER(len=200):: file_output


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CORPO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Lettura del file di input e del file contenente la mesh del problema
  CALL time2D_input(file_output)
  
  !Calcolo dei nodi e pesi della formula di quadratura di Gauss-Legendre
  CALL nodi_pesi_gauss

  !Costruzione dei blocchi della matrice di Toeplitz e dei corrispettivi
  !blocchi del termine noto
  CALL time2D_toeplitz_RHS(file_output)
 
  !Risoluzione del sistema e scrittura dell'output
  !CALL time2D_solver(file_output)
 
  !Calcolo della dimensione della matrice con allocazione di matrix e Vu
  !CALL DimSistema
  
  !Post-processing
  !CALL time2D_postpro(file_output)
  
END PROGRAM MAIN



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


