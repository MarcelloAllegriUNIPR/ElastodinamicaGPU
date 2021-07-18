SUBROUTINE time2D_postpro(file_output)
! Post processing

  USE variable_2Dgeneral

  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABILI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  CHARACTER(len=200),INTENT(IN):: file_output
  
  !Variabili locali
  INTEGER(kind=4):: i_time, ih, ik, iT, ind_el, kxspazio, kyspazio, tau, indice_termine_noto, iTT, ind_sol
  INTEGER(kind=4):: AllocateStatus
  
  REAL(kind=8):: x_min, x_max, y_min, y_max, trasl, xest, yest, deltax
  REAL(kind=8):: estremox1, estremoy1, estremox2, estremoy2, sumtau
  
  REAL(kind=8),DIMENSION(Nt)::sumcsi
  REAL(kind=8),DIMENSION(DimVu)::suu
  
  REAL(kind=8),DIMENSION(:),ALLOCATABLE:: griglia_x, griglia_y
  
  REAL(kind=8),EXTERNAL::Unota, DDOT
  
  CHARACTER(len=8):: iT_st, indice_termine_noto_st
  CHARACTER(len=8):: fmt2 !Descrittore di formato fmt
  
  CHARACTER(len=200)::  file_mesh_pp, file_sol !,ext, dir
  CHARACTER(len=200)::  file_path_mesh_pp, file_sol_pp, file_path_sol_pp, file_path_soluzione
  
  !!!!!!!!!!!!!!!!!!!!!!!!!! CORPO della SUBROUTINE !!!!!!!!!!!!!!!!!!

  !Calcolo delle massime e minime ascisse ed ordinate del dominio
  x_min=MINVAL(list_nodes(:)%coordinates(1),number_nodes)
  x_max=MAXVAL(list_nodes(:)%coordinates(1),number_nodes)
  y_min=MINVAL(list_nodes(:)%coordinates(2),number_nodes)
  y_max=MAXVAL(list_nodes(:)%coordinates(2),number_nodes)
  
  !Valore necessario per traslare assime e minime ascisse ed ordinate del dominio
  trasl=velC_P*T_fin+deltat

  !Ascisse ed ordinate degli estremi del dominio in cui effettuare il post processing
  !estremox1=x_min-trasl
  !estremox2=x_max+trasl
  !estremoy1=y_min-trasl
  !estremoy2=y_max+trasl
  
  estremox1=0.d0
  estremox2=10.d0
  estremoy1=0.d0
  estremoy2=0.d0

  deltax=velC_P*deltat !/2.d0 !passo di discretizzazione griglia esterna
  kxspazio=int((estremox2-estremox1)/deltax)
  kyspazio=int((estremoy2-estremoy1)/deltax)
  kyspazio=kxspazio
  
  !Allocazione della matrice contenente la griglia di punti per il post processing
  ALLOCATE(griglia_x(kxspazio+1),STAT=AllocateStatus)
  ALLOCATE(griglia_y(kyspazio+1),STAT=AllocateStatus)
  IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
  !Creazione della griglia
  DO ik=1,kxspazio+1
    griglia_x(ik)=(ik-1)*deltax+estremox1
  END DO
  DO ih=1,kyspazio+1
    !griglia_y(ih)=estremoy2-(ih-1)*deltax
	!griglia_y(ih)=estremoy1+(ih-1)*deltax
	 !griglia_y(ih)=dsin(3.d0*datan(1.d0))*griglia_x(ih)+0.5d0*dsin(datan(1.d0))
	 griglia_y(ih)=dsin(4.d0*datan(1.d0)/2.d0)*griglia_x(ih)
  END DO  
  DO ik=1,kxspazio+1
    !griglia_x(ik)=dcos(3.d0*datan(1.d0))*griglia_x(ik)-0.5d0+0.5d0*cos(datan(1.d0))
	griglia_x(ik)=dcos(4.d0*datan(1.d0)/2.d0)*griglia_x(ih)
  END DO
  
  !Formato per la conversione da integer a character della variabile i_time
  !(l'intero viene trasformato in una stringa di lunghezza 5 con 0 a sx)
  !fmt = '(I5.5)'
  !Formato per la conversione da integer a character della variabile indice_termine_noto
  !(l'intero viene scritto senza zeri (quindi 1 o 2))
  fmt2= '(I1)'
  !Nome della directory in cui si salva il file contenente il post-processing
  !dir='./../tests/output'
  !Estensione del file
  !ext='.txt' 
  
  !Nome del file
  file_mesh_pp=TRIM(file_output) // '_mesh_pp_piatto' 
  !Costruzione del path per il file di output in cui scrivere la mesh del post processing
  CALL filepath(file_path_mesh_pp,dir,file_mesh_pp,ext)
  !Apertura dei file di output
  OPEN(1,FILE=TRIM(file_path_mesh_pp)) 
  !WRITE(1,*) (kxspazio+1)*(kyspazio+1)
  !DO ik=1,kxspazio+1
    DO ih=1,kyspazio+1
        WRITE(1,*) griglia_x(ih), griglia_y(ih), ih !, ih
    END DO
  !END DO
  CLOSE(1)
  
  file_sol=TRIM(file_output) // '_soluzione'
  CALL filepath(file_path_soluzione,dir,file_sol,ext)
  !Allocazione del vettore soluzione
  ALLOCATE(sol(2*dimVu*Nt),STAT=AllocateStatus)!!!!!!!!!!!!!!!!!---MIA MODIFICA
  IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
  
  OPEN(10,FILE=file_path_soluzione,STATUS='OLD')
  !Ci posizioniamo all'inizio del file di input
  REWIND(10)! fa partire l'istruzione READ dalla prima riga del file

  !Lettura dell'intestazione del file mesh
  DO ind_sol=1,2*dimVu*Nt
     READ(10,*) sol(ind_sol)
  END DO
  CLOSE(10)
  
  write(*,*) 'discretizzazione spaziale', DimVu
  DO indice_termine_noto=1,2
  write(*,*) 'indice termine noto pp', indice_termine_noto
  DO iTT=1,8
    iT=iTT*20
  !DO it=1,Nt
    write(*,*) iT  
    !Conversione di iT da integer a character
    WRITE (iT_st,fmt) iT
	WRITE (indice_termine_noto_st,fmt2) indice_termine_noto
    !Nome del file
    file_sol_pp=TRIM(file_output) // '_sol_pp_piatto_' // TRIM(indice_termine_noto_st) // '_' // TRIM(iT_st)
    !Costruzione del path per il file di output in cui scrivere il post processing
    !all'istante iT
    CALL filepath(file_path_sol_pp,dir,file_sol_pp,ext)
    !Apertura dei file di output
    OPEN((indice_termine_noto-1)*Nt+iT,FILE=TRIM(file_path_sol_pp))!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !Valutazione del post processing all'istante iT
	DO ik=1,kxspazio+1
	ih=ik
	!DO ih=1,kyspazio+1
	    !if(griglia_y(ih).lt.1.d0)then
		     !sumtau=0.d0
		!(((0.d0.lt.(x_max+1.d-12)).and.(0.d0.gt.(x_min-1.d-12))).and.&!!originale !else
		if(((griglia_x(ik).lt.(x_max+1.d-12)).and.(griglia_x(ik).gt.(x_min-1.d-12))).and.&
		((griglia_y(ih).lt.(y_max+1.d-12)).and.(griglia_y(ih).gt.(y_min-1.d-12)))) then		
		    call dist(griglia_x(ik),griglia_y(ih),ind_el)				
		    if (ind_el.ge.1) then
				sumtau=Unota(list_elements(ind_el)%ind_BD,griglia_x(ik),griglia_y(ih),iT*deltat,indice_termine_noto)
				!write(*,*) 'faccio con la Unota'
				!write(*,*) 'risultato', sumtau
				!pause
		    else    
		        sumtau=0.d0
			    DO tau=0,Nt-1		
			        call MakeSolest(suu,griglia_x(ik),griglia_y(ih),iT*deltat,velC_S,velC_P,tau,indice_termine_noto,1)
                    !REF: http://www.netlib.org/lapack/explore-html/de/da4/group__double__blas__level1_ga75066c4825cb6ff1c8ec4403ef8c843a.html#ga75066c4825cb6ff1c8ec4403ef8c843a
					sumtau=sumtau+DDOT(DimVu,suu,1,sol(tau*2*DimVu+1:(tau+1)*2*DimVu-DimVu),1)
                    call MakeSolest(suu,griglia_x(ik),griglia_y(ih),iT*deltat,velC_S,velC_P,tau,indice_termine_noto,2)
                    !REF: http://www.netlib.org/lapack/explore-html/de/da4/group__double__blas__level1_ga75066c4825cb6ff1c8ec4403ef8c843a.html#ga75066c4825cb6ff1c8ec4403ef8c843a
					sumtau=sumtau+DDOT(DimVu,suu,1,sol(tau*2*DimVu+DimVu+1:(tau+1)*2*DimVu),1)		
        		END DO
				!write(*,*) 'faccio con proc 1'
				!write(*,*) 'risultato', sumtau
		    endif
		else
    		sumtau=0.d0
			DO tau=0,Nt-1			
			        call MakeSolest(suu,griglia_x(ik),griglia_y(ih),iT*deltat,velC_S,velC_P,tau,indice_termine_noto,1)
                    !REF: http://www.netlib.org/lapack/explore-html/de/da4/group__double__blas__level1_ga75066c4825cb6ff1c8ec4403ef8c843a.html#ga75066c4825cb6ff1c8ec4403ef8c843a
					sumtau=sumtau+DDOT(DimVu,suu,1,sol(tau*2*DimVu+1:(tau+1)*2*DimVu-DimVu),1)
                    call MakeSolest(suu,griglia_x(ik),griglia_y(ih),iT*deltat,velC_S,velC_P,tau,indice_termine_noto,2)
                    !REF: http://www.netlib.org/lapack/explore-html/de/da4/group__double__blas__level1_ga75066c4825cb6ff1c8ec4403ef8c843a.html#ga75066c4825cb6ff1c8ec4403ef8c843a
					sumtau=sumtau+DDOT(DimVu,suu,1,sol(tau*2*DimVu+DimVu+1:(tau+1)*2*DimVu),1)		
            END DO
			!write(*,*) 'faccio con proc 2'
			!write(*,*) 'risultato', sumtau
		endif
        WRITE((indice_termine_noto-1)*Nt+iT,*) sumtau, griglia_x(ik), griglia_y(ih)  !sumtau,  griglia_y(ih)
    END DO
    !END DO
	CLOSE((indice_termine_noto-1)*Nt+iT)
	!pause

  END DO 
  END DO
  
  RETURN
  END SUBROUTINE time2D_postpro

