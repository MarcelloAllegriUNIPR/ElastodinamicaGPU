SUBROUTINE time2D_read_input_file(file_path,file_mesh,file_out)

  USE variable_2Dgeneral

  IMPLICIT NONE

  !Input Parameters
  CHARACTER(len=200),INTENT(IN):: file_path

  !Output Parameters
  CHARACTER(len=200),INTENT(OUT):: file_mesh, file_out
  
  !Variabili Locali della subroutine
  CHARACTER(len=200):: rec

  !!!!!!!!!!!!!!!!!!!!!!!!!! CORPO della SUBROUTINE !!!!!!!!!!!!!!!!!!

  !Apertura del file di input
  OPEN(10,FILE=file_path,STATUS='old')! old vuol dire che il file già esiste (new invece crea un nuovo file)

  !Lettura del tipo di problema: frequency o time-domain
  READ(10,*) rec
  READ(10,*) rec

  !Lettura della dimensione dello spazio in cui è ambientato
  !il problema
  READ(10,*) rec
  READ(10,*) rec

  !Lettura del nome del file contenente la geometria della mesh
  READ(10,*) rec
  READ(10,*) file_mesh

  !Lettura del nome dei file di output
  READ(10,*) rec
  READ(10,*) file_out

  !Lettura del tempo finale
  READ(10,*) rec
  READ(10,*) T_fin
  
  !Lettura del numero di intervalli temporali
  READ(10,*) rec
  READ(10,*) Nt
  write(*,*) 'Nt', Nt
  
  !Calcolo del passo temporale
  deltat=T_fin/Nt
  write(*,*) 'deltat', deltat
  pause 
  
  !Lettura del grado delle funzioni di forma per approssimare
  !l'incognita u
  READ(10,*) rec
  READ(10,*) grado_u
  
  !Lettura del grado delle funzioni di forma per approssimare
  !l'incognita q
  READ(10,*) rec
  READ(10,*) grado_q

  !Lettura della densità di massa rho
  READ(10,*) rec
  READ(10,*) rho
  write(*,*) 'valore di rho', rho
  
  !Lettura del modulo di shear mu
  READ(10,*) rec
  READ(10,*) mu

  !Lettura del rapporto di Poisson nu
  READ(10,*) rec
  READ(10,*) nu

  !Calcolo del parametro di Lamé lambda
  !lambda=2.d0*mu*nu/(1.d0-2.d0*nu)
  !lambda=-mu    
  lambda=(2.d0)**2*rho-2.d0*mu
  write(*,*) 'valore di lambda', lambda
  pause
  
  !Calcolo della velocità delle onde P
  velC_P=DSQRT((lambda+2.d0*mu)/rho)
  write(*,*) 'cp', velC_P
  pause
  
  !Calcolo della velocità delle onde P
  velC_S=DSQRT(mu/rho)
  write(*,*) 'cs', velC_S
  pause
  
  !Chiusura del file di input
  CLOSE(10)
  
  RETURN
  
END SUBROUTINE time2D_read_input_file
