      Double precision function Vuextra_sub(l_m_tilde,l_m,e_m_tilde,e_m,tempo1,tempo2,indice_i,indice_j)
!
!     INTEGRALE DELLE FUNZIONI DI FORMA li,lj SUGLI ELEMENTI i,j
!     ELEMENTI SEPARATI

  USE variable_2Dgeneral
  use cudafor
  use OMP_LIB
  use VuExtraGpu
  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABILI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  INTEGER(kind=4),INTENT(IN):: e_m_tilde, e_m, l_m_tilde, l_m, indice_i, indice_j
  REAL(kind=8),INTENT(IN):: tempo1, tempo2
  
  !Variabili locali
  INTEGER(kind=4) :: ind_gauss, Ngauss, AllocateStatus, iplog, iqlog, ki, kj, Ngaussi, ii, jj, ierrSync
   
  REAL(kind=8) :: delta_x, alfa, beta, p2, A1, B1, alfa_j1, beta_j1, p2a, s, ds, xtrasl, xinttrasl, r2_1, serv, sm1, cs, cp, app
  
  REAL(kind=8),DIMENSION(:),ALLOCATABLE :: x, w, xint, wint, ti, vi, ww1
  
  REAL(kind=8),EXTERNAL :: fi1, dfi1, fiU!,Vuextra_sub_gpu
 
  INTEGER(kind=4):: flag_extra
  REAL(kind=8):: estremo_m, estremo_m_tilde, &
                 coeff_delta_kronecker, CA, CB, CC, CE, CF, CD, &
				 CBCFCECC, CFCACCCD, CBCCCECF, CACCCFCD, CACBCDCE, CBCDCECA, deltaquartiS, deltaquartiP, &
				 Vuextra_sub_P, Vuextra_sub_S, deltaquartiS_bis, deltaquartiP_bis, deltaquartiS_bbis, deltaquartiP_bbis, radice_incriminata
  INTEGER(kind=4), DIMENSION(2,2) :: delta_kronecker
  REAL(kind=8), DIMENSION(6) :: xx
  REAL(kind=8), DIMENSION(2) :: r, punto_m_1, punto_m_2, punto_m_tilde_1, punto_m_tilde_2

  DOUBLE PRECISION START1, END, value,  result

  type(dim3) :: dimGrid, dimBlock
  integer :: sizeInBytes, istat
  double precision, device :: alfa_d, beta_d, finalp2
  integer, device :: iplog_d, iqlog_d, blockNumber
  type(cudaEvent) :: start, stop
  real :: time
  !!!!!!!!!!!!!!!!!!!!!!!!!! CORPO della SUBROUTINE !!!!!!!!!!!!!!!!!!
  istat = cudaEventCreate(start)
  istat = cudaEventCreate(stop)
  
  cs=velC_S
  cp=velC_P
  
  punto_m_1=list_nodes(e_m)%coordinates(1:2)
  punto_m_2=list_nodes(e_m+1)%coordinates(1:2)
  punto_m_tilde_1=list_nodes(e_m_tilde)%coordinates(1:2)
  punto_m_tilde_2=list_nodes(e_m_tilde+1)%coordinates(1:2)
  
  delta_kronecker(1,1)=1
  delta_kronecker(1,2)=0
  delta_kronecker(2,1)=0
  delta_kronecker(2,2)=1
  coeff_delta_kronecker=delta_kronecker(indice_i,indice_j)/2.d0
  
Vuextra_sub=0.d0

delta_x=tempo1-tempo2
IF (delta_x.le.0.d0) RETURN

! IF (e_m_tilde.eq.81.and.e_m.eq.78.and.delta_kronecker(indice_i,indice_j).eq.0)then!!!!!!!!!!!!!!!!!!!!!!!!
  ! write(*,*) 'blocco', indice_i, indice_j!!!!!!!!!!!!!!!!!!
  ! write(*,*) 'tempo di chiamata', delta_x!!!!!!!!!!!!!!!!!!!
  ! pause!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Vuextra_sub_S=0.d0
  Vuextra_sub_P=0.d0
				      
  ind_gauss=6
  Ngauss=2**(ind_gauss-1) !32 nodi di Gauss

  ALLOCATE(x(Ngauss),STAT=AllocateStatus)
  IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
  ALLOCATE(w(Ngauss),STAT=AllocateStatus)
  IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
  ALLOCATE(xint(Ngauss),STAT=AllocateStatus)
  IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
  ALLOCATE(wint(Ngauss),STAT=AllocateStatus)
  IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    
  x=gauss(ind_gauss)%nodiquad
  w=gauss(ind_gauss)%pesiquad
  xint=gauss(ind_gauss)%nodiquad
  wint=gauss(ind_gauss)%pesiquad

  iplog=1
  iqlog=1

  estremo_m=sqrt((punto_m_2(1)-punto_m_1(1))**2+(punto_m_2(2)-punto_m_1(2))**2)
  estremo_m_tilde=sqrt((punto_m_tilde_2(1)-punto_m_tilde_1(1))**2+(punto_m_tilde_2(2)-punto_m_tilde_1(2))**2)
  
  CA=punto_m_tilde_1(1)-punto_m_1(1)
  CB=(punto_m_tilde_2(1)-punto_m_tilde_1(1))/estremo_m_tilde
  CC=(punto_m_2(1)-punto_m_1(1))/estremo_m
  CD=punto_m_tilde_1(2)-punto_m_1(2)
  CE=(punto_m_tilde_2(2)-punto_m_tilde_1(2))/estremo_m_tilde
  CF=(punto_m_2(2)-punto_m_1(2))/estremo_m

  CBCFCECC=CB*CF-CE*CC
  CFCACCCD=CF*CA-CC*CD
  CBCDCECA=CB*CD-CE*CA
  CBCCCECF=CB*CC+CE*CF
  CACCCFCD=CA*CC+CF*CD
  CACBCDCE=CA*CB+CD*CE 

  deltaquartiS=(cs*delta_x)**2-CFCACCCD**2
  deltaquartiP=(cp*delta_x)**2-CFCACCCD**2
  
  deltaquartiS_bis=(cs*delta_x)**2-CBCDCECA**2
  deltaquartiP_bis=(cp*delta_x)**2-CBCDCECA**2
  
  deltaquartiS_bbis=-(CBCDCECA)**2-estremo_m**2*(CBCFCECC)**2+2*estremo_m*CBCDCECA*CBCFCECC+(cs*delta_x)**2
  deltaquartiP_bbis=-(CBCDCECA)**2-estremo_m**2*(CBCFCECC)**2+2*estremo_m*CBCDCECA*CBCFCECC+(cp*delta_x)**2
     
  if((dabs(CBCFCECC).le.1.d-15).and.(dabs(CFCACCCD).le.1.d-15)) then
     flag_extra=1    !	elementi extra ma allineati
  elseif((dabs(CBCFCECC).le.1.d-15).and.(dabs(CFCACCCD).gt.1.d-15)) then
     flag_extra=2    !	paralleli 
  else
     flag_extra=3    !   non allineati non paralleli
  endif

  !wint, w, xint, x, cp, cs,grado_q
  if(useGpu .eq. 1) then
  call setInstanceCommonData(delta_x, indice_i, indice_j, CA, CB, CC, CD, CE, CF, &
                     coeff_delta_kronecker, flag_extra, estremo_m, l_m_tilde,l_m, estremo_m_tilde, &
                     CBCFCECC, CFCACCCD, CBCCCECF, CACCCFCD, CACBCDCE, CBCDCECA)
  endif
  !*************************************
  !             INTEGRAZIONE SU ES
  !*************************************
  IF((CBCFCECC.eq.0.d0).and.(deltaquartiS.lt.0.d0))then
     Vuextra_sub_S=0.d0
  else
     select case(flag_extra)
         case(1) 
	         if(CBCCCECF.lt.0.d0)then
                 xx(1)=-CACBCDCE+estremo_m*(CBCCCECF)-cs*delta_x
                 xx(4)=-CACBCDCE+cs*delta_x
                 xx(2)=dmin1(-CACBCDCE+estremo_m*(CBCCCECF)+cs*delta_x,-CACBCDCE-cs*delta_x)
		         xx(3)=dmax1(-CACBCDCE+estremo_m*(CBCCCECF)+cs*delta_x,-CACBCDCE-cs*delta_x)
			     xx(5)=xx(4)
			     xx(6)=xx(4)
		     else
		 	     xx(1)=-CACBCDCE-cs*delta_x
                 xx(4)=-CACBCDCE+estremo_m*(CBCCCECF)+cs*delta_x
                 xx(2)=dmin1(-CACBCDCE+estremo_m*(CBCCCECF)-cs*delta_x,-CACBCDCE+cs*delta_x)
		         xx(3)=dmax1(-CACBCDCE+estremo_m*(CBCCCECF)-cs*delta_x,-CACBCDCE+cs*delta_x)       
			     xx(5)=xx(4)
			     xx(6)=xx(4)
             endif
	     case(2)
		     if(CBCCCECF.lt.0.d0)then
                 xx(1)=-CACBCDCE+estremo_m*(CBCCCECF)-sqrt(deltaquartiS_bis)
                 xx(4)=-CACBCDCE+sqrt(deltaquartiS_bis)
				 xx(2)=dmin1(-CACBCDCE+estremo_m*(CBCCCECF)+sqrt(deltaquartiS_bis),-CACBCDCE-sqrt(deltaquartiS_bis))
				 xx(3)=dmax1(-CACBCDCE+estremo_m*(CBCCCECF)+sqrt(deltaquartiS_bis),-CACBCDCE-sqrt(deltaquartiS_bis))
				 xx(5)=xx(4)
				 xx(6)=xx(4)
			 else
                 xx(1)=-CACBCDCE-sqrt(deltaquartiS_bis)
	             xx(4)=-CACBCDCE+estremo_m*(CBCCCECF)+sqrt(deltaquartiS_bis)
				 xx(2)=dmin1(-CACBCDCE+estremo_m*(CBCCCECF)-sqrt(deltaquartiS_bis),-CACBCDCE+sqrt(deltaquartiS_bis))
				 xx(3)=dmax1(-CACBCDCE+estremo_m*(CBCCCECF)-sqrt(deltaquartiS_bis),-CACBCDCE+sqrt(deltaquartiS_bis))
				 xx(5)=xx(4)
				 xx(6)=xx(4)				 
             endif
         case(3)
		     xx(1)=dmin1((CFCACCCD-cs*delta_x)/(-CBCFCECC),(CFCACCCD+cs*delta_x)/(-CBCFCECC))
			 xx(6)=dmax1((CFCACCCD-cs*delta_x)/(-CBCFCECC),(CFCACCCD+cs*delta_x)/(-CBCFCECC))
			 if(deltaquartiS_bis.ge.0.d0)then
			     xx(2)=-(CACBCDCE)-sqrt(deltaquartiS_bis)
				 xx(3)=-(CACBCDCE)+sqrt(deltaquartiS_bis)
		     else		 
		     	 xx(2)=xx(1)
				 xx(3)=xx(1)
			 endif
			 if(deltaquartiS_bbis.ge.0.d0)then
			     xx(4)=-(CACBCDCE)+estremo_m*CBCCCECF-sqrt(deltaquartiS_bbis)
				 xx(5)=-(CACBCDCE)+estremo_m*CBCCCECF+sqrt(deltaquartiS_bbis)
		     else		 
		     	 xx(4)=xx(6)
				 xx(5)=xx(6)
			 endif
			 DO ii=2,6 !!!!!!!!!!!!!!!!ordinamento dei punti con insertion sort
                 app=xx(ii)       
                 jj=ii-1
                 do while (jj.ge.1)
                     if (xx(jj).gt.app) then
                         xx(jj+1)=xx(jj)
                         jj=jj-1 
                         xx(jj+1)=app
                     else 
                         exit
                     endif
	             enddo
             ENDDO
	 end select
	 do ii=1,6!!!!!!!!!!!!!!!!!taglio fuori i punti  minori di zero e maggiori di estremo_m_tilde
	     if(xx(ii).lt.0.d0)then
		     xx(ii)=0.d0
		 endif
         if(xx(ii).gt.estremo_m_tilde)then
		     xx(ii)=estremo_m_tilde
		 endif
	 enddo
	 DO ii=1,5!!!!!!!!!!!!!!!!!!!!!!!!integrazione nucleo su ES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	     !if(xx(ii+1)-xx(ii).gt.1.d-15)then
         if(xx(ii+1)-xx(ii).gt.1.d-14)then		 
			 alfa=(xx(ii+1)-xx(ii))/2.d0
	         beta=(xx(ii+1)+xx(ii))/2.d0
	
             p2 = 0.d0
	 
	         iplog=1
			 iqlog=1
             if(curva_piu_meno(beta,flag_extra,1,cs).le.estremo_m)then
			     iqlog=2
             elseif((curva_piu_meno(x(ii),flag_extra,1,cs).eq.estremo_m).or.(curva_piu_meno(x(ii+1),flag_extra,1,cs).eq.estremo_m))then
                 iqlog=2         
			 endif			

             if(curva_piu_meno(beta,flag_extra,-1,cs).ge.0.d0)then
			     iplog=2
             elseif((curva_piu_meno(x(ii),flag_extra,-1,cs).eq.0.d0).or.(curva_piu_meno(x(ii+1),flag_extra,-1,cs).eq.0.d0))then
                 iplog=2
			 endif

            if(useGpu .eq. 1) then
                alfa_d = alfa
                beta_d = beta
                iplog_d = iplog
                iqlog_d = iqlog

                dimGrid = dim3(Ngauss,1,1)
                dimBlock = dim3(1,Ngauss,1)             
                sizeInBytes = sizeof(app)*Ngauss

                istat = cudaEventRecord(start,0)
                call performCalcES<<<dimGrid,dimBlock,sizeInBytes>>>(alfa_d, beta_d,iplog_d, iqlog_d)
                
                dimGrid = dim3(1,1,1)
                dimBlock = dim3(Ngauss,1,1)
                finalp2 = 0.d0
                call finalSum<<<dimGrid,dimBlock>>>(finalp2)
                istat = cudaEventRecord(stop,0)
                istat = cudaDeviceSynchronize()
                istat = cudaEventElapsedTime(time, start, stop)
                p2 = finalp2
                gputime = gputime + time/(1.0e3)

                ierrSync = cudaGetLastError()
                if (ierrSync /= cudaSuccess) then
                    write(*,*) 'Sync kernel error:', cudaGetErrorString(ierrSync)
                endif

                !iterationCounter = iterationCounter + 1
                !print *, "iterationCounter", iterationCounter
                !print *, p2
                !pause
            else
                p2 = 0.d0
                START1 = omp_get_wtime() 
    	     	DO ki=1,Ngauss	
                 	xtrasl=alfa*x(ki)+beta
				    A1=dmax1(0.d0,curva_piu_meno(xtrasl,flag_extra,-1,cs))
           	    	B1=dmin1(estremo_m,curva_piu_meno(xtrasl,flag_extra,1,cs))
				 	alfa_j1=(B1-A1)/2.d0
	             	beta_j1=(B1+A1)/2.d0
				 	p2a = 0.d0
				 	IF ((B1-A1).gt.10.d-14) THEN	                 
                    	DO kj=1,Ngauss
	                     	xinttrasl=(xint(kj)+1.d0)*0.5d0
	                     	s=fi1(iplog,iqlog,xinttrasl)
	                     	ds=dfi1(iplog,iqlog,xinttrasl)
	                     	serv=alfa_j1*(2.d0*s-1.d0)+beta_j1
                         	r2_1=(CA+CB*xtrasl-CC*serv)**2+(CD+CE*xtrasl-CF*serv)**2
						 	r(1)=CA+CB*xtrasl-CC*serv
						 	r(2)=CD+CE*xtrasl-CF*serv
	                     	p2a = p2a - wint(kj)*ds*fiU(l_m,serv,estremo_m,grado_q)*(r(indice_i)*r(indice_j)/(r2_1**2)-coeff_delta_kronecker/r2_1)*(delta_x/cs)*sqrt(dabs((cs*delta_x)**2-r2_1))							
	                     	IF (delta_kronecker(indice_i,indice_j).eq.1.d0) THEN
	                         	p2a= p2a+wint(kj)*ds*fiU(l_m,serv,estremo_m,grado_q)*coeff_delta_kronecker*(1/cs**2)*(dlog(cs*delta_x+sqrt(dabs((cs*delta_x)**2-r2_1)))-dlog(sqrt(r2_1)))
							ENDIF							
                    	END DO
        	     	ENDIF
	             	p2 = p2 +p2a*alfa_j1*w(ki)*fiU(l_m_tilde,xtrasl,estremo_m_tilde,grado_q)                    
                    
                    !value = p2a*alfa_j1*w(ki)*fiU(l_m_tilde,xtrasl,estremo_m_tilde,grado_q)
                    !p2 = p2+value
                    !print *,"cpu", ki, value
                END DO
                cputime =  cputime + omp_get_wtime() - START1
                !print *,"cpu",p2
                !if(abs(value-p2) .gt. 1.d-10) then
                !    print *, "ES cpu", p2, "gpu", value
                !    pause
                !endif
            endif            
            !pause
            Vuextra_sub_S=Vuextra_sub_S+p2*alfa             
		else
		    Vuextra_sub_S=Vuextra_sub_S+0.d0
        endif		 
	 enddo	 
  endif
  
  ! write(*,*) 'integrale cs', Vuextra_sub_S!!!!!!!!!!!!!!!!!
  ! pause!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  
  !*************************************
  !             INTEGRAZIONE SU EP
  !*************************************
  IF((CBCFCECC.eq.0.d0).and.(deltaquartiP.lt.0.d0))then
     Vuextra_sub_P=0.d0
  else
     select case(flag_extra)
         case(1) 
	         if(CBCCCECF.lt.0.d0)then
                 xx(1)=-CACBCDCE+estremo_m*(CBCCCECF)-cp*delta_x
                 xx(4)=-CACBCDCE+cp*delta_x
                 xx(2)=dmin1(-CACBCDCE+estremo_m*(CBCCCECF)+cp*delta_x,-CACBCDCE-cp*delta_x)
		         xx(3)=dmax1(-CACBCDCE+estremo_m*(CBCCCECF)+cp*delta_x,-CACBCDCE-cp*delta_x)
			     xx(5)=xx(4)
			     xx(6)=xx(4)
		     else
		 	     xx(1)=-CACBCDCE-cp*delta_x
                 xx(4)=-CACBCDCE+estremo_m*(CBCCCECF)+cp*delta_x
                 xx(2)=dmin1(-CACBCDCE+estremo_m*(CBCCCECF)-cp*delta_x,-CACBCDCE+cp*delta_x)
		         xx(3)=dmax1(-CACBCDCE+estremo_m*(CBCCCECF)-cp*delta_x,-CACBCDCE+cp*delta_x)       
			     xx(5)=xx(4)
			     xx(6)=xx(4)
             endif
	     case(2)
		     if(CBCCCECF.lt.0.d0)then
                 xx(1)=-CACBCDCE+estremo_m*(CBCCCECF)-sqrt(deltaquartiP_bis)
                 xx(4)=-CACBCDCE+sqrt(deltaquartiP_bis)
				 xx(2)=dmin1(-CACBCDCE+estremo_m*(CBCCCECF)+sqrt(deltaquartiP_bis),-CACBCDCE-sqrt(deltaquartiP_bis))
				 xx(3)=dmax1(-CACBCDCE+estremo_m*(CBCCCECF)+sqrt(deltaquartiP_bis),-CACBCDCE-sqrt(deltaquartiP_bis))
				 xx(5)=xx(4)
				 xx(6)=xx(4)
			 else
                 xx(1)=-CACBCDCE-sqrt(deltaquartiP_bis)
	             xx(4)=-CACBCDCE+estremo_m*(CBCCCECF)+sqrt(deltaquartiP_bis)
				 xx(2)=dmin1(-CACBCDCE+estremo_m*(CBCCCECF)-sqrt(deltaquartiP_bis),-CACBCDCE+sqrt(deltaquartiP_bis))
				 xx(3)=dmax1(-CACBCDCE+estremo_m*(CBCCCECF)-sqrt(deltaquartiP_bis),-CACBCDCE+sqrt(deltaquartiP_bis))
				 xx(5)=xx(4)
				 xx(6)=xx(4)				 
             endif
         case(3)
		     xx(1)=dmin1((CFCACCCD-cp*delta_x)/(-CBCFCECC),(CFCACCCD+cp*delta_x)/(-CBCFCECC))
			 xx(6)=dmax1((CFCACCCD-cp*delta_x)/(-CBCFCECC),(CFCACCCD+cp*delta_x)/(-CBCFCECC))
			 if(deltaquartiP_bis.ge.0.d0)then
			     xx(2)=-(CACBCDCE)-sqrt(deltaquartiP_bis)
				 xx(3)=-(CACBCDCE)+sqrt(deltaquartiP_bis)
		     else		 
		     	 xx(2)=xx(1)
				 xx(3)=xx(1)
			 endif
			 if(deltaquartiP_bbis.ge.0.d0)then
			     xx(4)=-(CACBCDCE)+estremo_m*CBCCCECF-sqrt(deltaquartiP_bbis)
				 xx(5)=-(CACBCDCE)+estremo_m*CBCCCECF+sqrt(deltaquartiP_bbis)
		     else		 
		     	 xx(4)=xx(6)
				 xx(5)=xx(6)
			 endif
			 DO ii=2,6 !!!!!!!!!!!!!!!!ordinamento dei punti con insertion sort
                 app=xx(ii)       
                 jj=ii-1
                 do while (jj.ge.1)
                     if (xx(jj).gt.app) then
                         xx(jj+1)=xx(jj)
                         jj=jj-1 
                         xx(jj+1)=app
                     else 
                         exit
                     endif
	             enddo
             ENDDO
	 end select
	 do ii=1,6!!!!!!!!!!!!!!!!!taglio fuori i punti  minori di zero e maggiori di estremo_m_tilde
	     if(xx(ii).lt.0.d0)then
		     xx(ii)=0.d0
		 endif
         if(xx(ii).gt.estremo_m_tilde)then
		     xx(ii)=estremo_m_tilde
		 endif
	 enddo
	 DO ii=1,5!!!!!!!!!!!!!!!!!!!!!!!!integrazione nucleo su EP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	     !if(xx(ii+1)-xx(ii).gt.1.d-15)then
		 if(xx(ii+1)-xx(ii).gt.1.d-14)then
			 alfa=(xx(ii+1)-xx(ii))/2.d0
	         beta=(xx(ii+1)+xx(ii))/2.d0
	
             p2 = 0.d0

	         iplog=1
			 iqlog=1
             if(curva_piu_meno(beta,flag_extra,1,cp).le.estremo_m)then
			     iqlog=2
             elseif((curva_piu_meno(x(ii),flag_extra,1,cp).eq.estremo_m).or.(curva_piu_meno(x(ii+1),flag_extra,1,cp).eq.estremo_m))then
                 iqlog=2         
			 endif			

             if(curva_piu_meno(beta,flag_extra,-1,cp).ge.0.d0)then
			     iplog=2
             elseif((curva_piu_meno(x(ii),flag_extra,-1,cp).eq.0.d0).or.(curva_piu_meno(x(ii+1),flag_extra,-1,cp).eq.0.d0))then
                 iplog=2
			 endif

            if(useGpu .eq. 1) then
                alfa_d = alfa
                beta_d = beta
                iplog_d = iplog
                iqlog_d = iqlog

                dimGrid = dim3(Ngauss,1,1)
                dimBlock = dim3(1,Ngauss,1)             
                sizeInBytes = sizeof(app)*Ngauss

                istat = cudaEventRecord(start,0)
                call performCalcEP<<<dimGrid,dimBlock,sizeInBytes>>>(alfa_d, beta_d,iplog_d, iqlog_d)
                
                dimGrid = dim3(1,1,1)
                dimBlock = dim3(Ngauss,1,1)
                finalp2 = 0.d0
                call finalSum<<<dimGrid,dimBlock>>>(finalp2)
                istat = cudaEventRecord(stop,0)
                istat = cudaDeviceSynchronize()
                istat = cudaEventElapsedTime(time, start, stop)
                p2 = finalp2
                gputime = gputime + time/(1.0e3)

                ierrSync = cudaGetLastError()
                if (ierrSync /= cudaSuccess) then
                    write(*,*) 'Sync kernel error:', cudaGetErrorString(ierrSync)
                endif

                !iterationCounter = iterationCounter + 1
                !print *, "iterationCounter", iterationCounter
                !print *, p2
                !pause
            else
                p2 = 0.d0
                START1 = omp_get_wtime() 
                DO ki=1,Ngauss	
                    xtrasl=alfa*x(ki)+beta
                    A1=dmax1(0.d0,curva_piu_meno(xtrasl,flag_extra,-1,cp))
                    B1=dmin1(estremo_m,curva_piu_meno(xtrasl,flag_extra,1,cp))	             
                    alfa_j1=(B1-A1)/2.d0
                    beta_j1=(B1+A1)/2.d0
                    p2a = 0.d0		                    
                    !IF ((B1-A1).gt.10.d-15) THEN
                    IF ((B1-A1).gt.10.d-14) THEN	
                        DO kj=1,Ngauss
                            !value = 0.d0
                            xinttrasl=(xint(kj)+1.d0)*0.5d0
                            s=fi1(iplog,iqlog,xinttrasl)
                            ds=dfi1(iplog,iqlog,xinttrasl)
                            serv=alfa_j1*(2.d0*s-1.d0)+beta_j1
                            r2_1=(CA+CB*xtrasl-CC*serv)**2+(CD+CE*xtrasl-CF*serv)**2
                            r(1)=CA+CB*xtrasl-CC*serv
                            r(2)=CD+CE*xtrasl-CF*serv
                            p2a = p2a+wint(kj)*ds*fiU(l_m,serv,estremo_m,grado_q)*(r(indice_i)*r(indice_j)/(r2_1**2)-coeff_delta_kronecker/r2_1)*(delta_x/cp)*sqrt(dabs((cp*delta_x)**2-r2_1))
	                        !p2a = p2a + wint(kj)*ds*fiU(l_m,serv,estremo_m,grado_q)*(r(indice_i)*r(indice_j)/(r2_1**2)-coeff_delta_kronecker/r2_1)*(delta_x/cp)*sqrt((cp*delta_x)**2-r2_1)
	                        
                            !if(ki .eq. 16) then
                            !    print *, kj, wint(kj)
                            !endif

                            IF (delta_kronecker(indice_i,indice_j).eq.1.d0) THEN
                                p2a= p2a+wint(kj)*ds*fiU(l_m,serv,estremo_m,grado_q)*coeff_delta_kronecker*(1/cp**2)*(dlog(cp*delta_x+sqrt(dabs((cp*delta_x)**2-r2_1)))-dlog(sqrt(r2_1)))
                                !p2a= p2a+wint(kj)*ds*fiU(l_m,serv,estremo_m,grado_q)*coeff_delta_kronecker*(1/cp**2)*(dlog(cp*delta_x+sqrt((cp*delta_x)**2-r2_1))-dlog(sqrt(r2_1)))
                            ENDIF	
                            !p2a = p2a + value
                            !if(ki .eq. 16) then
                            !    print *, kj, value
                            !endif                            
                        END DO
                        !if(ki .eq. 16) then
                        !    print *, ki, p2a
                        !endif
                    ENDIF
                    p2=p2+p2a*alfa_j1*w(ki)*fiU(l_m_tilde,xtrasl,estremo_m_tilde,grado_q)
                    !value = p2a*alfa_j1*w(ki)*fiU(l_m_tilde,xtrasl,estremo_m_tilde,grado_q)
                    !p2 = p2 + value                    
                    !print *, "cpu p2:",ki, value
                    !print *,"cpu", ki, p2a                    				 
                END DO
                cputime =  cputime + omp_get_wtime() - START1
                
                !iterationCounter = iterationCounter + 1
                !print *, "iterationCounter", iterationCounter
                !print *,p2
                !pause
                !if(abs(value-p2) .gt. 1.d-10) then
                !    print *, "EP cpu", p2, "gpu", value
                !    pause
                !endif               
            endif            
             Vuextra_sub_P=Vuextra_sub_P+p2*alfa		 
		 else
		     Vuextra_sub_P=Vuextra_sub_P+0.d0
         endif		 
	 enddo	 
  endif
  Vuextra_sub=Vuextra_sub_P+Vuextra_sub_S
  !print *, Vuextra_sub_P,Vuextra_sub_S
  !pause
  ! write(*,*) 'integrale cp', Vuextra_sub_P!!!!!!!!!!!!!!!!!
  ! write(*,*) 'integrale', Vuextra_sub!!!!!!!!!!!!!!!!!
  ! pause!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  
  ! else!!!!!!!!!!!!!!!!!!!!!!!!
     ! Vuextra_sub=0.d0!!!!!!!!!!!!!!!!!!!!
  ! endif!!!!!!!!!!!!!!!!!!!!!!!!!
  
RETURN

   CONTAINS 
   
   FUNCTION curva_piu_meno(x_linea,tipo_allineamento,segno,vel) !!! per definire le rette dei domini in caso di elementi allineati (o eventualmete limiti delle ordinate)
      REAL(kind=8) :: curva_piu_meno, x_linea, vel
      INTEGER(kind=4) :: tipo_allineamento, segno
	  select case (tipo_allineamento)
         case (1)!!!!!allineati
             curva_piu_meno=CBCCCECF*x_linea+CACCCFCD+sign(1,segno)*vel*delta_x
         case (2)!!!!!paralleli
             curva_piu_meno=CBCCCECF*x_linea+CACCCFCD+sign(1,segno)*sqrt(dabs((vel*delta_x)**2-(CFCACCCD)**2))
         case (3)!!!!!generici
             curva_piu_meno=CBCCCECF*x_linea+CACCCFCD+sign(1,segno)*sqrt(dabs(-(CBCFCECC*x_linea)**2+2*CBCFCECC*(-CFCACCCD)*x_linea-(CFCACCCD)**2+vel**2*delta_x**2))
      end select
   END FUNCTION curva_piu_meno
END