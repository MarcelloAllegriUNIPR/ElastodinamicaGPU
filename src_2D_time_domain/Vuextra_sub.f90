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
 
  INTEGER(kind=4):: flag_extra_vuextra
  REAL(kind=8):: estremo_m_Vuextra, estremo_m_tilde_Vuextra, &
                 coeff_delta_kronecker_vuextra, CA_vuextra, CB_vuextra, CC_vuextra, CE_vuextra, CF_vuextra, CD_vuextra, &
				 CBCFCECC_vuextra, CFCACCCD_vuextra, CBCCCECF_vuextra, CACCCFCD_vuextra, CACBCDCE_vuextra, CBCDCECA_vuextra, deltaquartiS, deltaquartiP, &
				 Vuextra_sub_P, Vuextra_sub_S, deltaquartiS_bis, deltaquartiP_bis, deltaquartiS_bbis, deltaquartiP_bbis, radice_incriminata
  INTEGER(kind=4), DIMENSION(2,2) :: delta_kronecker
  REAL(kind=8), DIMENSION(6) :: xx
  REAL(kind=8), DIMENSION(2) :: r_vuextra, punto_m_1_vuextra, punto_m_2_vuextra, punto_m_tilde_1_vuextra, punto_m_tilde_2_vuextra

  DOUBLE PRECISION START1, END, value

  integer :: istat
  type(cudaEvent) :: start, stop  
  real :: time
!   double precision, device :: alfa_d, beta_d, result !CalculationResults(10),
!   integer, device :: iplog_d, iqlog_d
  !!!!!!!!!!!!!!!!!!!!!!!!!! CORPO della SUBROUTINE !!!!!!!!!!!!!!!!!!
  istat = cudaEventCreate(start)
  istat = cudaEventCreate(stop)
  
  cs=velC_S
  cp=velC_P
  
  punto_m_1_vuextra=list_nodes(e_m)%coordinates(1:2)
  punto_m_2_vuextra=list_nodes(e_m+1)%coordinates(1:2)
  punto_m_tilde_1_vuextra=list_nodes(e_m_tilde)%coordinates(1:2)
  punto_m_tilde_2_vuextra=list_nodes(e_m_tilde+1)%coordinates(1:2)
  
  delta_kronecker(1,1)=1
  delta_kronecker(1,2)=0
  delta_kronecker(2,1)=0
  delta_kronecker(2,2)=1
  coeff_delta_kronecker_vuextra=delta_kronecker(indice_i,indice_j)/2.d0
  
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

  estremo_m_Vuextra=sqrt((punto_m_2_vuextra(1)-punto_m_1_vuextra(1))**2+(punto_m_2_vuextra(2)-punto_m_1_vuextra(2))**2)
  estremo_m_tilde_Vuextra=sqrt((punto_m_tilde_2_vuextra(1)-punto_m_tilde_1_vuextra(1))**2+(punto_m_tilde_2_vuextra(2)-punto_m_tilde_1_vuextra(2))**2)
  
  CA_vuextra=punto_m_tilde_1_vuextra(1)-punto_m_1_vuextra(1)
  CB_vuextra=(punto_m_tilde_2_vuextra(1)-punto_m_tilde_1_vuextra(1))/estremo_m_tilde_Vuextra
  CC_vuextra=(punto_m_2_vuextra(1)-punto_m_1_vuextra(1))/estremo_m_Vuextra
  CD_vuextra=punto_m_tilde_1_vuextra(2)-punto_m_1_vuextra(2)
  CE_vuextra=(punto_m_tilde_2_vuextra(2)-punto_m_tilde_1_vuextra(2))/estremo_m_tilde_Vuextra
  CF_vuextra=(punto_m_2_vuextra(2)-punto_m_1_vuextra(2))/estremo_m_Vuextra

  CBCFCECC_vuextra=CB_vuextra*CF_vuextra-CE_vuextra*CC_vuextra
  CFCACCCD_vuextra=CF_vuextra*CA_vuextra-CC_vuextra*CD_vuextra
  CBCDCECA_vuextra=CB_vuextra*CD_vuextra-CE_vuextra*CA_vuextra
  CBCCCECF_vuextra=CB_vuextra*CC_vuextra+CE_vuextra*CF_vuextra
  CACCCFCD_vuextra=CA_vuextra*CC_vuextra+CF_vuextra*CD_vuextra
  CACBCDCE_vuextra=CA_vuextra*CB_vuextra+CD_vuextra*CE_vuextra 

  deltaquartiS=(cs*delta_x)**2-CFCACCCD_vuextra**2
  deltaquartiP=(cp*delta_x)**2-CFCACCCD_vuextra**2
  
  deltaquartiS_bis=(cs*delta_x)**2-CBCDCECA_vuextra**2
  deltaquartiP_bis=(cp*delta_x)**2-CBCDCECA_vuextra**2
  
  deltaquartiS_bbis=-(CBCDCECA_vuextra)**2-estremo_m_Vuextra**2*(CBCFCECC_vuextra)**2+2*estremo_m_Vuextra*CBCDCECA_vuextra*CBCFCECC_vuextra+(cs*delta_x)**2
  deltaquartiP_bbis=-(CBCDCECA_vuextra)**2-estremo_m_Vuextra**2*(CBCFCECC_vuextra)**2+2*estremo_m_Vuextra*CBCDCECA_vuextra*CBCFCECC_vuextra+(cp*delta_x)**2
     
  if((dabs(CBCFCECC_vuextra).le.1.d-15).and.(dabs(CFCACCCD_vuextra).le.1.d-15)) then
     flag_extra_vuextra=1    !	elementi extra ma allineati
  elseif((dabs(CBCFCECC_vuextra).le.1.d-15).and.(dabs(CFCACCCD_vuextra).gt.1.d-15)) then
     flag_extra_vuextra=2    !	paralleli 
  else
     flag_extra_vuextra=3    !   non allineati non paralleli
  endif

  !wint, w, xint, x, cp, cs,grado_q
  if(useGpu .eq. 1) then
    ! call setInstanceCommonData(delta_x, indice_i, indice_j, CA_vuextra, CB_vuextra, CC_vuextra, CD_vuextra, CE_vuextra, CF_vuextra, &
    !                  coeff_delta_kronecker_vuextra, flag_extra_vuextra, estremo_m_Vuextra, l_m_tilde,l_m, estremo_m_tilde_Vuextra, &
    !                  CBCFCECC_vuextra, CFCACCCD_vuextra, CBCCCECF_vuextra, CACCCFCD_vuextra, CACBCDCE_vuextra, CBCDCECA_vuextra)
    !CalculationResults = 0.d0  
     !                result = 0.d0
endif
  !*************************************
  !             INTEGRAZIONE SU ES
  !*************************************
  IF((CBCFCECC_vuextra.eq.0.d0).and.(deltaquartiS.lt.0.d0))then
     Vuextra_sub_S=0.d0
  else
     select case(flag_extra_vuextra)
         case(1) 
	         if(CBCCCECF_vuextra.lt.0.d0)then
                 xx(1)=-CACBCDCE_vuextra+estremo_m_Vuextra*(CBCCCECF_vuextra)-cs*delta_x
                 xx(4)=-CACBCDCE_vuextra+cs*delta_x
                 xx(2)=dmin1(-CACBCDCE_vuextra+estremo_m_Vuextra*(CBCCCECF_vuextra)+cs*delta_x,-CACBCDCE_vuextra-cs*delta_x)
		         xx(3)=dmax1(-CACBCDCE_vuextra+estremo_m_Vuextra*(CBCCCECF_vuextra)+cs*delta_x,-CACBCDCE_vuextra-cs*delta_x)
			     xx(5)=xx(4)
			     xx(6)=xx(4)
		     else
		 	     xx(1)=-CACBCDCE_vuextra-cs*delta_x
                 xx(4)=-CACBCDCE_vuextra+estremo_m_Vuextra*(CBCCCECF_vuextra)+cs*delta_x
                 xx(2)=dmin1(-CACBCDCE_vuextra+estremo_m_Vuextra*(CBCCCECF_vuextra)-cs*delta_x,-CACBCDCE_vuextra+cs*delta_x)
		         xx(3)=dmax1(-CACBCDCE_vuextra+estremo_m_Vuextra*(CBCCCECF_vuextra)-cs*delta_x,-CACBCDCE_vuextra+cs*delta_x)       
			     xx(5)=xx(4)
			     xx(6)=xx(4)
             endif
	     case(2)
		     if(CBCCCECF_vuextra.lt.0.d0)then
                 xx(1)=-CACBCDCE_vuextra+estremo_m_Vuextra*(CBCCCECF_vuextra)-sqrt(deltaquartiS_bis)
                 xx(4)=-CACBCDCE_vuextra+sqrt(deltaquartiS_bis)
				 xx(2)=dmin1(-CACBCDCE_vuextra+estremo_m_Vuextra*(CBCCCECF_vuextra)+sqrt(deltaquartiS_bis),-CACBCDCE_vuextra-sqrt(deltaquartiS_bis))
				 xx(3)=dmax1(-CACBCDCE_vuextra+estremo_m_Vuextra*(CBCCCECF_vuextra)+sqrt(deltaquartiS_bis),-CACBCDCE_vuextra-sqrt(deltaquartiS_bis))
				 xx(5)=xx(4)
				 xx(6)=xx(4)
			 else
                 xx(1)=-CACBCDCE_vuextra-sqrt(deltaquartiS_bis)
	             xx(4)=-CACBCDCE_vuextra+estremo_m_Vuextra*(CBCCCECF_vuextra)+sqrt(deltaquartiS_bis)
				 xx(2)=dmin1(-CACBCDCE_vuextra+estremo_m_Vuextra*(CBCCCECF_vuextra)-sqrt(deltaquartiS_bis),-CACBCDCE_vuextra+sqrt(deltaquartiS_bis))
				 xx(3)=dmax1(-CACBCDCE_vuextra+estremo_m_Vuextra*(CBCCCECF_vuextra)-sqrt(deltaquartiS_bis),-CACBCDCE_vuextra+sqrt(deltaquartiS_bis))
				 xx(5)=xx(4)
				 xx(6)=xx(4)				 
             endif
         case(3)
		     xx(1)=dmin1((CFCACCCD_vuextra-cs*delta_x)/(-CBCFCECC_vuextra),(CFCACCCD_vuextra+cs*delta_x)/(-CBCFCECC_vuextra))
			 xx(6)=dmax1((CFCACCCD_vuextra-cs*delta_x)/(-CBCFCECC_vuextra),(CFCACCCD_vuextra+cs*delta_x)/(-CBCFCECC_vuextra))
			 if(deltaquartiS_bis.ge.0.d0)then
			     xx(2)=-(CACBCDCE_vuextra)-sqrt(deltaquartiS_bis)
				 xx(3)=-(CACBCDCE_vuextra)+sqrt(deltaquartiS_bis)
		     else		 
		     	 xx(2)=xx(1)
				 xx(3)=xx(1)
			 endif
			 if(deltaquartiS_bbis.ge.0.d0)then
			     xx(4)=-(CACBCDCE_vuextra)+estremo_m_Vuextra*CBCCCECF_vuextra-sqrt(deltaquartiS_bbis)
				 xx(5)=-(CACBCDCE_vuextra)+estremo_m_Vuextra*CBCCCECF_vuextra+sqrt(deltaquartiS_bbis)
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
	 do ii=1,6!!!!!!!!!!!!!!!!!taglio fuori i punti  minori di zero e maggiori di estremo_m_tilde_Vuextra
	     if(xx(ii).lt.0.d0)then
		     xx(ii)=0.d0
		 endif
         if(xx(ii).gt.estremo_m_tilde_Vuextra)then
		     xx(ii)=estremo_m_tilde_Vuextra
		 endif
	 enddo
     if(useGpu .eq. 0) then
        DO ii=1,5!!!!!!!!!!!!!!!!!!!!!!!!integrazione nucleo su ES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !if(xx(ii+1)-xx(ii).gt.1.d-15)then            
            if(xx(ii+1)-xx(ii).gt.1.d-14)then		 
                alfa=(xx(ii+1)-xx(ii))/2.d0
                beta=(xx(ii+1)+xx(ii))/2.d0
        
                p2 = 0.d0
        
                iplog=1
                iqlog=1
                if(curva_piu_meno(beta,flag_extra_vuextra,1,cs).le.estremo_m_Vuextra)then
                    iqlog=2
                elseif((curva_piu_meno(x(ii),flag_extra_vuextra,1,cs).eq.estremo_m_Vuextra).or.(curva_piu_meno(x(ii+1),flag_extra_vuextra,1,cs).eq.estremo_m_Vuextra))then
                    iqlog=2         
                endif			

                if(curva_piu_meno(beta,flag_extra_vuextra,-1,cs).ge.0.d0)then
                    iplog=2
                elseif((curva_piu_meno(x(ii),flag_extra_vuextra,-1,cs).eq.0.d0).or.(curva_piu_meno(x(ii+1),flag_extra_vuextra,-1,cs).eq.0.d0))then
                    iplog=2
                endif
                
                p2 = 0.d0
                START1 = omp_get_wtime()
                DO ki=1,Ngauss	
                    xtrasl=alfa*x(ki)+beta
                    A1=dmax1(0.d0,curva_piu_meno(xtrasl,flag_extra_vuextra,-1,cs))
                    B1=dmin1(estremo_m_Vuextra,curva_piu_meno(xtrasl,flag_extra_vuextra,1,cs))
                    alfa_j1=(B1-A1)/2.d0
                    beta_j1=(B1+A1)/2.d0
                    p2a = 0.d0                    
                    IF ((B1-A1).gt.10.d-14) THEN	                 
                        DO kj=1,Ngauss     
                            ! value = 0.d0                       
                            xinttrasl=(xint(kj)+1.d0)*0.5d0
                            s=fi1(iplog,iqlog,xinttrasl)
                            ds=dfi1(iplog,iqlog,xinttrasl)
                            serv=alfa_j1*(2.d0*s-1.d0)+beta_j1
                            r2_1=(CA_vuextra+CB_vuextra*xtrasl-CC_vuextra*serv)**2+(CD_vuextra+CE_vuextra*xtrasl-CF_vuextra*serv)**2
                            r_vuextra(1)=CA_vuextra+CB_vuextra*xtrasl-CC_vuextra*serv
                            r_vuextra(2)=CD_vuextra+CE_vuextra*xtrasl-CF_vuextra*serv
                            p2a = p2a - wint(kj)*ds*fiU(l_m,serv,estremo_m_Vuextra,grado_q)*(r_vuextra(indice_i)*r_vuextra(indice_j)/(r2_1**2)-coeff_delta_kronecker_vuextra/r2_1)*(delta_x/cs)*sqrt(dabs((cs*delta_x)**2-r2_1))							                            
                            !if(ki .eq. 16) value = - wint(kj)*ds*fiU(l_m,serv,estremo_m_Vuextra,grado_q)*(r_vuextra(indice_i)*r_vuextra(indice_j)/(r2_1**2)-coeff_delta_kronecker_vuextra/r2_1)*(delta_x/cs)*sqrt(dabs((cs*delta_x)**2-r2_1))
                            !value = -wint(kj)*ds*fiU(l_m,serv,estremo_m_Vuextra,grado_q)*(r_vuextra(indice_i)*r_vuextra(indice_j)/(r2_1**2)-coeff_delta_kronecker_vuextra/r2_1)*(delta_x/cs)*sqrt(dabs((cs*delta_x)**2-r2_1))							                            
                            IF (delta_kronecker(indice_i,indice_j).eq.1.d0) THEN
                                p2a= p2a+wint(kj)*ds*fiU(l_m,serv,estremo_m_Vuextra,grado_q)*coeff_delta_kronecker_vuextra*(1/cs**2)*(dlog(cs*delta_x+sqrt(dabs((cs*delta_x)**2-r2_1)))-dlog(sqrt(r2_1)))                                
                                !value= value+wint(kj)*ds*fiU(l_m,serv,estremo_m_Vuextra,grado_q)*coeff_delta_kronecker_vuextra*(1/cs**2)*(dlog(cs*delta_x+sqrt(dabs((cs*delta_x)**2-r2_1)))-dlog(sqrt(r2_1)))                                
                                ! if(ki .eq. 16) then
                                !     value= value+wint(kj)*ds*fiU(l_m,serv,estremo_m_Vuextra,grado_q)*coeff_delta_kronecker_vuextra*(1/cs**2)*(dlog(cs*delta_x+sqrt(dabs((cs*delta_x)**2-r2_1)))-dlog(sqrt(r2_1)))
                                !     print *,value, ki,kj
                                ! endif
                            ENDIF                            					
                        END DO
                    ENDIF
                    p2 = p2 + p2a*alfa_j1*w(ki)*fiU(l_m_tilde,xtrasl,estremo_m_tilde_Vuextra,grado_q)
                    ! print *,value, ki
                    !  +value                                 
                END DO
                cputime =  cputime + omp_get_wtime() - START1               
                !pause
                ! print *,p2*alfa, '267'
                Vuextra_sub_S=Vuextra_sub_S+p2*alfa             
            else
                Vuextra_sub_S=Vuextra_sub_S+0.d0
            endif		 
        enddo     
    else
        DO ii=1,5!!!!!!!!!!!!!!!!!!!!!!!!integrazione nucleo su ES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !if(xx(ii+1)-xx(ii).gt.1.d-15)then
            if(xx(ii+1)-xx(ii).gt.1.d-14)then		 
                ! alfa_d=(xx(ii+1)-xx(ii))/2.d0
                ! beta_d=(xx(ii+1)+xx(ii))/2.d0

                ! iplog_d=1
                ! iqlog_d=1
                ! if(curva_piu_meno(beta,flag_extra_vuextra,1,cs).le.estremo_m_Vuextra)then
                !     iqlog_d=2
                ! elseif((curva_piu_meno(x(ii),flag_extra_vuextra,1,cs).eq.estremo_m_Vuextra).or.(curva_piu_meno(x(ii+1),flag_extra_vuextra,1,cs).eq.estremo_m_Vuextra))then
                !     iqlog_d=2         
                ! endif			

                ! if(curva_piu_meno(beta,flag_extra_vuextra,-1,cs).ge.0.d0)then
                !     iplog_d=2
                ! elseif((curva_piu_meno(x(ii),flag_extra_vuextra,-1,cs).eq.0.d0).or.(curva_piu_meno(x(ii+1),flag_extra_vuextra,-1,cs).eq.0.d0))then
                !     iplog_d=2
                ! endif
                
                !istat = cudaEventRecord(start,0)
                
                !print *,"call preCalculationES<<<dimGrid,dimBlockPreCalculation>>>(alfa_d,beta_d,iplog_d, iqlog_d)"
                !call preCalculationES<<<dimGrid,dimBlockPreCalculation>>>(alfa_d,beta_d,iplog_d, iqlog_d)
                !print *,"call performCalcES<<<dimGrid,dimBlockCalculation>>>(alfa_d,result)"
                !call performCalcES<<<dimGrid,dimBlockCalculation,sizeof(app)*NGaussDimension*NGaussDimension>>>(alfa_d,result)                
                !Vuextra_sub = Vuextra_sub + result
                !pause

                !istat = cudaEventRecord(stop,0)
                !istat = cudaDeviceSynchronize()
                !istat = cudaEventElapsedTime(time, start, stop)                
                !gputime = gputime + time/(1.0e3)

                ! ierrSync = cudaGetLastError()
                ! if (ierrSync /= cudaSuccess) then
                !     write(*,*) 'Sync kernel error:', cudaGetErrorString(ierrSync)
                ! endif
            endif
        enddo  
    endif
  endif
  
  ! write(*,*) 'integrale cs', Vuextra_sub_S!!!!!!!!!!!!!!!!!
  ! pause!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  
  !*************************************
  !             INTEGRAZIONE SU EP
  !*************************************
  IF((CBCFCECC_vuextra.eq.0.d0).and.(deltaquartiP.lt.0.d0))then
     Vuextra_sub_P=0.d0
  else
     select case(flag_extra_vuextra)
         case(1) 
	         if(CBCCCECF_vuextra.lt.0.d0)then
                 xx(1)=-CACBCDCE_vuextra+estremo_m_Vuextra*(CBCCCECF_vuextra)-cp*delta_x
                 xx(4)=-CACBCDCE_vuextra+cp*delta_x
                 xx(2)=dmin1(-CACBCDCE_vuextra+estremo_m_Vuextra*(CBCCCECF_vuextra)+cp*delta_x,-CACBCDCE_vuextra-cp*delta_x)
		         xx(3)=dmax1(-CACBCDCE_vuextra+estremo_m_Vuextra*(CBCCCECF_vuextra)+cp*delta_x,-CACBCDCE_vuextra-cp*delta_x)
			     xx(5)=xx(4)
			     xx(6)=xx(4)
		     else
		 	     xx(1)=-CACBCDCE_vuextra-cp*delta_x
                 xx(4)=-CACBCDCE_vuextra+estremo_m_Vuextra*(CBCCCECF_vuextra)+cp*delta_x
                 xx(2)=dmin1(-CACBCDCE_vuextra+estremo_m_Vuextra*(CBCCCECF_vuextra)-cp*delta_x,-CACBCDCE_vuextra+cp*delta_x)
		         xx(3)=dmax1(-CACBCDCE_vuextra+estremo_m_Vuextra*(CBCCCECF_vuextra)-cp*delta_x,-CACBCDCE_vuextra+cp*delta_x)       
			     xx(5)=xx(4)
			     xx(6)=xx(4)
             endif
	     case(2)
		     if(CBCCCECF_vuextra.lt.0.d0)then
                 xx(1)=-CACBCDCE_vuextra+estremo_m_Vuextra*(CBCCCECF_vuextra)-sqrt(deltaquartiP_bis)
                 xx(4)=-CACBCDCE_vuextra+sqrt(deltaquartiP_bis)
				 xx(2)=dmin1(-CACBCDCE_vuextra+estremo_m_Vuextra*(CBCCCECF_vuextra)+sqrt(deltaquartiP_bis),-CACBCDCE_vuextra-sqrt(deltaquartiP_bis))
				 xx(3)=dmax1(-CACBCDCE_vuextra+estremo_m_Vuextra*(CBCCCECF_vuextra)+sqrt(deltaquartiP_bis),-CACBCDCE_vuextra-sqrt(deltaquartiP_bis))
				 xx(5)=xx(4)
				 xx(6)=xx(4)
			 else
                 xx(1)=-CACBCDCE_vuextra-sqrt(deltaquartiP_bis)
	             xx(4)=-CACBCDCE_vuextra+estremo_m_Vuextra*(CBCCCECF_vuextra)+sqrt(deltaquartiP_bis)
				 xx(2)=dmin1(-CACBCDCE_vuextra+estremo_m_Vuextra*(CBCCCECF_vuextra)-sqrt(deltaquartiP_bis),-CACBCDCE_vuextra+sqrt(deltaquartiP_bis))
				 xx(3)=dmax1(-CACBCDCE_vuextra+estremo_m_Vuextra*(CBCCCECF_vuextra)-sqrt(deltaquartiP_bis),-CACBCDCE_vuextra+sqrt(deltaquartiP_bis))
				 xx(5)=xx(4)
				 xx(6)=xx(4)				 
             endif
         case(3)
		     xx(1)=dmin1((CFCACCCD_vuextra-cp*delta_x)/(-CBCFCECC_vuextra),(CFCACCCD_vuextra+cp*delta_x)/(-CBCFCECC_vuextra))
			 xx(6)=dmax1((CFCACCCD_vuextra-cp*delta_x)/(-CBCFCECC_vuextra),(CFCACCCD_vuextra+cp*delta_x)/(-CBCFCECC_vuextra))
			 if(deltaquartiP_bis.ge.0.d0)then
			     xx(2)=-(CACBCDCE_vuextra)-sqrt(deltaquartiP_bis)
				 xx(3)=-(CACBCDCE_vuextra)+sqrt(deltaquartiP_bis)
		     else		 
		     	 xx(2)=xx(1)
				 xx(3)=xx(1)
			 endif
			 if(deltaquartiP_bbis.ge.0.d0)then
			     xx(4)=-(CACBCDCE_vuextra)+estremo_m_Vuextra*CBCCCECF_vuextra-sqrt(deltaquartiP_bbis)
				 xx(5)=-(CACBCDCE_vuextra)+estremo_m_Vuextra*CBCCCECF_vuextra+sqrt(deltaquartiP_bbis)
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
	 do ii=1,6!!!!!!!!!!!!!!!!!taglio fuori i punti  minori di zero e maggiori di estremo_m_tilde_Vuextra
	     if(xx(ii).lt.0.d0)then
		     xx(ii)=0.d0
		 endif
         if(xx(ii).gt.estremo_m_tilde_Vuextra)then
		     xx(ii)=estremo_m_tilde_Vuextra
		 endif
	 enddo
    if(useGpu .eq. 0) then
        DO ii=1,5!!!!!!!!!!!!!!!!!!!!!!!!integrazione nucleo su EP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !if(xx(ii+1)-xx(ii).gt.1.d-15)then
            if(xx(ii+1)-xx(ii).gt.1.d-14)then
                alfa=(xx(ii+1)-xx(ii))/2.d0
                beta=(xx(ii+1)+xx(ii))/2.d0
        
                p2 = 0.d0

                iplog=1
                iqlog=1
                if(curva_piu_meno(beta,flag_extra_vuextra,1,cp).le.estremo_m_Vuextra)then
                    iqlog=2
                elseif((curva_piu_meno(x(ii),flag_extra_vuextra,1,cp).eq.estremo_m_Vuextra).or.(curva_piu_meno(x(ii+1),flag_extra_vuextra,1,cp).eq.estremo_m_Vuextra))then
                    iqlog=2         
                endif			

                if(curva_piu_meno(beta,flag_extra_vuextra,-1,cp).ge.0.d0)then
                    iplog=2
                elseif((curva_piu_meno(x(ii),flag_extra_vuextra,-1,cp).eq.0.d0).or.(curva_piu_meno(x(ii+1),flag_extra_vuextra,-1,cp).eq.0.d0))then
                    iplog=2
                endif

                
                
                    p2 = 0.d0
                    START1 = omp_get_wtime() 
                    DO ki=1,Ngauss	
                        xtrasl=alfa*x(ki)+beta
                        A1=dmax1(0.d0,curva_piu_meno(xtrasl,flag_extra_vuextra,-1,cp))
                        B1=dmin1(estremo_m_Vuextra,curva_piu_meno(xtrasl,flag_extra_vuextra,1,cp))	             
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
                                r2_1=(CA_vuextra+CB_vuextra*xtrasl-CC_vuextra*serv)**2+(CD_vuextra+CE_vuextra*xtrasl-CF_vuextra*serv)**2
                                r_vuextra(1)=CA_vuextra+CB_vuextra*xtrasl-CC_vuextra*serv
                                r_vuextra(2)=CD_vuextra+CE_vuextra*xtrasl-CF_vuextra*serv
                                !p2a = p2a+wint(kj)*ds*fiU(l_m,serv,estremo_m_Vuextra,grado_q)*(r_vuextra(indice_i)*r_vuextra(indice_j)/(r2_1**2)-coeff_delta_kronecker_vuextra/r2_1)*(delta_x/cp)*sqrt(dabs((cp*delta_x)**2-r2_1))
                                !value = +wint(kj)*ds*fiU(l_m,serv,estremo_m_Vuextra,grado_q)*(r_vuextra(indice_i)*r_vuextra(indice_j)/(r2_1**2)-coeff_delta_kronecker_vuextra/r2_1)*(delta_x/cp)*sqrt(dabs((cp*delta_x)**2-r2_1))
                                !p2a = p2a + wint(kj)*ds*fiU(l_m,serv,estremo_m_Vuextra,grado_q)*(r_vuextra(indice_i)*r_vuextra(indice_j)/(r2_1**2)-coeff_delta_kronecker_vuextra/r2_1)*(delta_x/cp)*sqrt((cp*delta_x)**2-r2_1)
                                
                                ! if(ki .eq. 16) then
                                !     print *, kj, wint(kj)
                                ! endif

                                IF (delta_kronecker(indice_i,indice_j).eq.1.d0) THEN
                                    p2a= p2a+wint(kj)*ds*fiU(l_m,serv,estremo_m_Vuextra,grado_q)*coeff_delta_kronecker_vuextra*(1/cp**2)*(dlog(cp*delta_x+sqrt(dabs((cp*delta_x)**2-r2_1)))-dlog(sqrt(r2_1)))
                                    !value= value+wint(kj)*ds*fiU(l_m,serv,estremo_m_Vuextra,grado_q)*coeff_delta_kronecker_vuextra*(1/cp**2)*(dlog(cp*delta_x+sqrt(dabs((cp*delta_x)**2-r2_1)))-dlog(sqrt(r2_1)))
                                    !p2a= p2a+wint(kj)*ds*fiU(l_m,serv,estremo_m_Vuextra,grado_q)*coeff_delta_kronecker_vuextra*(1/cp**2)*(dlog(cp*delta_x+sqrt((cp*delta_x)**2-r2_1))-dlog(sqrt(r2_1)))
                                ENDIF	
                                ! p2a = p2a + value
                                ! if(ki .eq. 16) then
                                !     print *, kj, value
                                ! endif                            
                            END DO
                            !if(ki .eq. 16) then
                            !    print *, ki, p2a
                            !endif
                        ENDIF
                        p2=p2+p2a*alfa_j1*w(ki)*fiU(l_m_tilde,xtrasl,estremo_m_tilde_Vuextra,grado_q)
                        !value = p2a*alfa_j1*w(ki)*fiU(l_m_tilde,xtrasl,estremo_m_tilde_Vuextra,grado_q)
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
                ! print *,p2*alfa,'479'
                Vuextra_sub_P=Vuextra_sub_P+p2*alfa		 
            else
                Vuextra_sub_P=Vuextra_sub_P+0.d0
            endif		 
        enddo	     
    else
        DO ii=1,5!!!!!!!!!!!!!!!!!!!!!!!!integrazione nucleo su EP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            
             if(xx(ii+1)-xx(ii).gt.1.d-14)then
            !     alfa_d=(xx(ii+1)-xx(ii))/2.d0
            !     beta_d=(xx(ii+1)+xx(ii))/2.d0
        
            !     iplog_d=1
            !     iqlog_d=1
            !     if(curva_piu_meno(beta,flag_extra_vuextra,1,cp).le.estremo_m_Vuextra)then
            !         iqlog_d=2
            !     elseif((curva_piu_meno(x(ii),flag_extra_vuextra,1,cp).eq.estremo_m_Vuextra).or.(curva_piu_meno(x(ii+1),flag_extra_vuextra,1,cp).eq.estremo_m_Vuextra))then
            !         iqlog_d=2         
            !     endif			

            !     if(curva_piu_meno(beta,flag_extra_vuextra,-1,cp).ge.0.d0)then
            !         iplog_d=2
            !     elseif((curva_piu_meno(x(ii),flag_extra_vuextra,-1,cp).eq.0.d0).or.(curva_piu_meno(x(ii+1),flag_extra_vuextra,-1,cp).eq.0.d0))then
            !         iplog_d=2
            !     endif
		 
                !istat = cudaEventRecord(start,0)                

                !print *, "call preCalculationEP<<<dimGrid,dimBlockPreCalculation>>>(alfa_d,beta_d,iplog_d, iqlog_d)"
                !call preCalculationEP<<<dimGrid,dimBlockPreCalculation>>>(alfa_d,beta_d,iplog_d, iqlog_d)
                !print *, "performCalcEP<<<dimGrid,dimBlockCalculation>>>(alfa_d,result)"                
                !call performCalcEP<<<dimGrid,dimBlockCalculation,sizeof(app)*NGaussDimension*NGaussDimension>>>(alfa_d,result)                
                !Vuextra_sub = Vuextra_sub + result

                !istat = cudaEventRecord(stop,0)
                !istat = cudaDeviceSynchronize()
                !istat = cudaEventElapsedTime(time, start, stop)                
                !gputime = gputime + time/(1.0e3)

                !ierrSync = cudaGetLastError()
                !if (ierrSync /= cudaSuccess) then
                !   write(*,*) 'Sync kernel error:', cudaGetErrorString(ierrSync)
                !endif            
            endif		 
        enddo	 
    endif
  endif
  !if(useGpu .eq. 0) then
      Vuextra_sub=Vuextra_sub_P+Vuextra_sub_S 
      !pause     
  !else
    
   ! do ii=1,10
        !value = CalculationResults(ii)
        !print *,value
        !Vuextra_sub = result
    !enddo    
  !endif
  
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
             curva_piu_meno=CBCCCECF_vuextra*x_linea+CACCCFCD_vuextra+sign(1,segno)*vel*delta_x
         case (2)!!!!!paralleli
             curva_piu_meno=CBCCCECF_vuextra*x_linea+CACCCFCD_vuextra+sign(1,segno)*sqrt(dabs((vel*delta_x)**2-(CFCACCCD_vuextra)**2))
         case (3)!!!!!generici
             curva_piu_meno=CBCCCECF_vuextra*x_linea+CACCCFCD_vuextra+sign(1,segno)*sqrt(dabs(-(CBCFCECC_vuextra*x_linea)**2+2*CBCFCECC_vuextra*(-CFCACCCD_vuextra)*x_linea-(CFCACCCD_vuextra)**2+vel**2*delta_x**2))
      end select
   END FUNCTION curva_piu_meno
END