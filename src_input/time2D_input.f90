SUBROUTINE time2D_input(file_output)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                    
! SUBROUTINE NAME: time2D_input
!
! GOAL: la subroutine legge il file di input del problema (./../tests/
!       input/input_time.txt), il quale ha la seguente struttura:
!
!       *\METODO
!       time
!       *\DIMENSIONE dello SPAZIO
!       2
!       *\NOME FILE MESH
!       nome_file_mesh 
!       *\NOME FILE di OUTPUT
!       nome_file_output 
!       *\TEMPO FINALE 
!       valore intero
!       *\NUMERO INTERVALLI TEMPORALI
!       valore intero
!       *\GRADO u
!       valore intero
!       *\GRADO q
!       valore intero
!       *\DENSITA’ di MASSA rho
!       valore double
!       *\SHEAR MODULO mu
!       valore double
!       *\RAPPORTO di POISSON nu 
!       valore double  
!
!       Successivamente la mesh del problema viene letta dal file
!       nome_file_mesh.mesh, a partire dalla quale vengono ricavate
!       tutte le informazioni geometriche necessarie per calcolare
!       gli integrali  
!                                                                    
! OUTPUT :  file_output, nome dei file di output
! 
! AUTORE: Luca Desiderio                                             
!                                                                    
! DATA:   29/XII/2017                           
!                                                                    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

  IMPLICIT NONE
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABILI !!!!!!!!!!!!!!!!!!!!!!!!!!

  !Output
  CHARACTER(len=200),INTENT(OUT):: file_output

  !Variabili locali della subroutine
  CHARACTER(len=200):: file_path, file_mesh, file, dir, ext
  
  !!!!!!!!!!!!!!!!!!!!!!!!!! CORPO della SUBROUTINE !!!!!!!!!!!!!!!!!!

  !***** STEP 1: file di input
  
  !Nome della directory in cui si trova il file di input
  dir='./../tests/input'
  !Nome del file di input
  file='input_time'
  !Estensione del file di input
  ext='.txt'

  !Costruzione del path per il file di input
  CALL filepath(file_path,dir,file,ext)

  !Verifica dell'esistenza del file di input
  CALL check_file(file_path)
  
  !Lettura del file di input
  CALL time2D_read_input_file(file_path,file_mesh,file_output)

  !***** STEP 2: file mesh

  !Nome della directory in cui si trova il file mesh
  dir='./../tests/mesh'
  !Estensione del file di input
  ext='.mesh'

  !Costruzione del path per il file mesh
  CALL filepath(file_path,dir,file_mesh,ext)

  !Verifica dell'esistenza del file mesh
  CALL check_file(file_path)

  !Lettura della mesh e rielaborazione dei dati
  CALL time2D_read_mesh(file_path)
  
  !***** STEP 3: file mesh (tipo ALBORGHETTI)
  
  !Nome della directory in cui si mettere il file mesh tipo ALBORGHETTI
  dir='./../tests/output'
  !Nome del file
  file=TRIM(file_output) // '_mesh_estesa'
  !Estensione del file di input
  ext='.txt'
  
  !Costruzione del path per il file di output
  CALL filepath(file_path,dir,file,ext)
  
  !Riscrittura file di dati
  CALL write_output_2D_descrizione_problema(file_path)

  RETURN
  
END SUBROUTINE time2D_input
