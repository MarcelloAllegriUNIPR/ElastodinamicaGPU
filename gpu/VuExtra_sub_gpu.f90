subroutine Vuextra_sub_gpu(CalculationResults,rowIndex)

  USE variable_2Dgeneral
  use variables
  use cudafor
  use OMP_LIB
  use VuExtraGpu
  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABILI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Double precision, intent(out), device :: CalculationResults(10,4)  
  integer, intent(in) :: rowIndex
  !Input
  !INTEGER(kind=4),INTENT(IN):: e_m_tilde, e_m, l_m_tilde, l_m, indice_i, indice_j
  !REAL(kind=8),INTENT(IN):: tempo1, tempo2
  
  REAL(kind=8) :: delta_x, alfa, beta, cs, cp, app,alfa_h(5)
  !Variabili locali
  INTEGER(kind=4) :: AllocateStatus, iplog, iqlog,position
  INTEGER(kind=4)::  ii, jj,ierrSync,ierrAsync !flag_extra
  REAL(kind=8):: deltaquartiS, deltaquartiP, Vuextra_sub_P, Vuextra_sub_S, deltaquartiS_bis, deltaquartiP_bis, deltaquartiS_bbis, deltaquartiP_bbis
  
  REAL(kind=8), DIMENSION(6) :: xx  
  DOUBLE PRECISION START1, END
  REAL(kind=8),EXTERNAL :: fi1, dfi1, fiU

  integer :: istat, ind_gauss = 6,ki
  type(cudaEvent) :: start, stop  
  double precision :: time
  double precision, device :: alfa_d(10), beta_d
  integer, device :: iplog_d, iqlog_d, rowIndex_d,position_d(10)
  !!!!!!!!!!!!!!!!!!!!!!!!!! CORPO della SUBROUTINE !!!!!!!!!!!!!!!!!!
  
delta_x=delta_x_VuExtra(rowIndex)
IF (delta_x.le.0.d0) RETURN  

    rowIndex_d = rowIndex
  deltaquartiS=(velC_S*delta_x)**2-CFCACCCD**2
  deltaquartiP=(velC_P*delta_x)**2-CFCACCCD**2
  
  deltaquartiS_bis=(velC_S*delta_x)**2-CBCDCECA**2
  deltaquartiP_bis=(velC_P*delta_x)**2-CBCDCECA**2
  
  deltaquartiS_bbis=-(CBCDCECA)**2-estremo_m**2*(CBCFCECC)**2+2*estremo_m*CBCDCECA*CBCFCECC+(velC_S*delta_x)**2
  deltaquartiP_bbis=-(CBCDCECA)**2-estremo_m**2*(CBCFCECC)**2+2*estremo_m*CBCDCECA*CBCFCECC+(velC_P*delta_x)**2 



  !*************************************
  !             INTEGRAZIONE SU ES
  !*************************************


  IF((CBCFCECC.eq.0.d0).and.(deltaquartiS.lt.0.d0))then
     Vuextra_sub_S=0.d0
  else
     select case(flag_extra)
         case(1) 
	         if(CBCCCECF.lt.0.d0)then
                 xx(1)=-CACBCDCE+estremo_m*(CBCCCECF)-velC_S*delta_x
                 xx(4)=-CACBCDCE+velC_S*delta_x
                 xx(2)=dmin1(-CACBCDCE+estremo_m*(CBCCCECF)+velC_S*delta_x,-CACBCDCE-velC_S*delta_x)
		         xx(3)=dmax1(-CACBCDCE+estremo_m*(CBCCCECF)+velC_S*delta_x,-CACBCDCE-velC_S*delta_x)
			     xx(5)=xx(4)
			     xx(6)=xx(4)
		     else
		 	     xx(1)=-CACBCDCE-velC_S*delta_x
                 xx(4)=-CACBCDCE+estremo_m*(CBCCCECF)+velC_S*delta_x
                 xx(2)=dmin1(-CACBCDCE+estremo_m*(CBCCCECF)-velC_S*delta_x,-CACBCDCE+velC_S*delta_x)
		         xx(3)=dmax1(-CACBCDCE+estremo_m*(CBCCCECF)-velC_S*delta_x,-CACBCDCE+velC_S*delta_x)       
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
		     xx(1)=dmin1((CFCACCCD-velC_S*delta_x)/(-CBCFCECC),(CFCACCCD+velC_S*delta_x)/(-CBCFCECC))
			 xx(6)=dmax1((CFCACCCD-velC_S*delta_x)/(-CBCFCECC),(CFCACCCD+velC_S*delta_x)/(-CBCFCECC))
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
            if(xx(ii+1)-xx(ii).gt.1.d-14)then		 
                alfa=(xx(ii+1)-xx(ii))/2.d0
                alfa_d(ii)=alfa
                beta=(xx(ii+1)+xx(ii))/2.d0

                !print *, "b"
                iplog=1
                iqlog=1
                if(curva_piu_meno(beta,flag_extra,1,velC_S).le.estremo_m)then
                    iqlog=2
                elseif((curva_piu_meno(x_VuExtra(ii),flag_extra,1,velC_S).eq.estremo_m).or.(curva_piu_meno(x_VuExtra(ii+1),flag_extra,1,velC_S).eq.estremo_m))then
                    iqlog=2         
                endif			

                if(curva_piu_meno(beta,flag_extra,-1,velC_S).ge.0.d0)then
                    iplog=2
                elseif((curva_piu_meno(x_VuExtra(ii),flag_extra,-1,velC_S).eq.0.d0).or.(curva_piu_meno(x_VuExtra(ii+1),flag_extra,-1,velC_S).eq.0.d0))then
                    iplog=2
                endif
                
                !call preCalculationES<<<dimGrid,dimBlockPreCalculation>>>(alfa_d,beta_d,iplog_d, iqlog_d)
                position = ((rowIndex-1)*10)+ii
                position_d(ii)=position
                
                !!$omp parallel do
                DO ki=1,32	
                    xtrasl(ki)=alfa*x_VuExtra(ki)+beta
                    A1(ki)=dmax1(0.d0,curva_piu_meno(xtrasl(ki),flag_extra,-1,velC_S))
                    B1(ki)=dmin1(estremo_m,curva_piu_meno(xtrasl(ki),flag_extra,1,velC_S))
                    alfa_j1(ki)=(B1(ki)-A1(ki))/2.d0
                    beta_j1(ki)=(B1(ki)+A1(ki))/2.d0                    
                    IF ((B1(ki)-A1(ki)).gt.10.d-14) THEN                                            
                        xinttrasl(ki)=(x_VuExtra(ki)+1.d0)*0.5d0
                        s(ki)=fi1(iplog,iqlog,xinttrasl(ki))
                        ds(ki)=dfi1(iplog,iqlog,xinttrasl(ki))
                    else            
                        xinttrasl(ki)=0.d0
                        s(ki)=0.d0
                        ds(ki)=0.d0                    
                    endif
                enddo
                !!$omp end parallel do

                !print *, rowIndex, ii, 'copy es'
                ! istat = cudaMemcpy2DAsync(xtrasl_d(1,position),32,xtrasl(1),32,32,1,1,stream(position))
                ! istat = cudaMemcpy2DAsync(A1_d(1,position),32,A1(1),32,32,1,1,stream(position))
                ! istat = cudaMemcpy2DAsync(B1_d(1,position),32,B1(1),32,32,1,1,stream(position))
                ! istat = cudaMemcpy2DAsync(alfa_j1_d(1,position),32,alfa_j1(1),32,32,1,1,stream(position))
                ! istat = cudaMemcpy2DAsync(beta_j1_d(1,position),32,beta_j1(1),32,32,1,1,stream(position))
                ! istat = cudaMemcpy2DAsync(xinttrasl_d(1,position),32,xtrasl(1),32,32,1,1,stream(position))
                ! istat = cudaMemcpy2DAsync(s_d(1,position),32,s(1),32,32,1,1,stream(position))
                ! istat = cudaMemcpy2DAsync(ds_d(1,position),32,ds(1),32,32,1,1,stream(position))
                
                ! call performCalcES<<<dimGrid,dimBlockCalculation,sharedMemDimension,stream(position)>>>(alfa_d(ii),position_d(ii),rowIndex_d)
                ! call sum_p2a<<<dimGrid,dimBlockCalculation,sharedMemDimension,stream(position)>>>(alfa_d(ii),position_d(ii), CalculationResults(ii,rowIndex))
                istat = cudaMemcpy2DAsync(xtrasl_d(1,position),32,xtrasl(1),32,32,1,1)
                istat = cudaMemcpy2DAsync(A1_d(1,position),32,A1(1),32,32,1,1)
                istat = cudaMemcpy2DAsync(B1_d(1,position),32,B1(1),32,32,1,1)
                istat = cudaMemcpy2DAsync(alfa_j1_d(1,position),32,alfa_j1(1),32,32,1,1)
                istat = cudaMemcpy2DAsync(beta_j1_d(1,position),32,beta_j1(1),32,32,1,1)
                istat = cudaMemcpy2DAsync(xinttrasl_d(1,position),32,xtrasl(1),32,32,1,1)
                istat = cudaMemcpy2DAsync(s_d(1,position),32,s(1),32,32,1,1)
                istat = cudaMemcpy2DAsync(ds_d(1,position),32,ds(1),32,32,1,1)
                
                call performCalcES<<<dimGrid,dimBlockCalculation,sharedMemDimension>>>(alfa_d(ii),position_d(ii),rowIndex_d)
                call sum_p2a<<<dimGrid,dimBlockCalculation,sharedMemDimension>>>(alfa_d(ii),position_d(ii), CalculationResults(ii,rowIndex))
            endif
        enddo
        ! do ii=1,5
        !     if(xx(ii+1)-xx(ii).gt.1.d-14)then
        !         position = ((rowIndex-1)*10)+ii         
                
        !         ! ierrSync = cudaGetLastError()
        !         ! ierrAsync = cudaDeviceSynchronize()
        !         ! if (ierrSync /= cudaSuccess) write(*,*) 'Sync kernel error:', cudaGetErrorString(ierrSync)
        !         ! if (ierrAsync /= cudaSuccess) write(*,*) 'Async kernel error:', cudaGetErrorString(ierrAsync)

        !         !pause            
        !         ! app = CalculationResults(ii,rowIndex)
        !         ! print *, app, ii, rowIndex
        !     endif
        ! enddo
        ! pause       
  endif
  
  ! write(*,*) 'integrale velC_S', Vuextra_sub_S!!!!!!!!!!!!!!!!!
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
                 xx(1)=-CACBCDCE+estremo_m*(CBCCCECF)-velC_P*delta_x
                 xx(4)=-CACBCDCE+velC_P*delta_x
                 xx(2)=dmin1(-CACBCDCE+estremo_m*(CBCCCECF)+velC_P*delta_x,-CACBCDCE-velC_P*delta_x)
		         xx(3)=dmax1(-CACBCDCE+estremo_m*(CBCCCECF)+velC_P*delta_x,-CACBCDCE-velC_P*delta_x)
			     xx(5)=xx(4)
			     xx(6)=xx(4)
		     else
		 	     xx(1)=-CACBCDCE-velC_P*delta_x
                 xx(4)=-CACBCDCE+estremo_m*(CBCCCECF)+velC_P*delta_x
                 xx(2)=dmin1(-CACBCDCE+estremo_m*(CBCCCECF)-velC_P*delta_x,-CACBCDCE+velC_P*delta_x)
		         xx(3)=dmax1(-CACBCDCE+estremo_m*(CBCCCECF)-velC_P*delta_x,-CACBCDCE+velC_P*delta_x)       
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
		     xx(1)=dmin1((CFCACCCD-velC_P*delta_x)/(-CBCFCECC),(CFCACCCD+velC_P*delta_x)/(-CBCFCECC))
			 xx(6)=dmax1((CFCACCCD-velC_P*delta_x)/(-CBCFCECC),(CFCACCCD+velC_P*delta_x)/(-CBCFCECC))
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
            if(xx(ii+1)-xx(ii).gt.1.d-14)then
                alfa=(xx(ii+1)-xx(ii))/2.d0
                alfa_d(ii+5)=alfa
                beta=(xx(ii+1)+xx(ii))/2.d0
        
                iplog=1
                iqlog=1
                if(curva_piu_meno(beta,flag_extra,1,velC_P).le.estremo_m)then
                    iqlog=2
                elseif((curva_piu_meno(x_VuExtra(ii),flag_extra,1,velC_P).eq.estremo_m).or.(curva_piu_meno(x_VuExtra(ii+1),flag_extra,1,velC_P).eq.estremo_m))then
                    iqlog=2         
                endif			

                if(curva_piu_meno(beta,flag_extra,-1,velC_P).ge.0.d0)then
                    iplog=2
                elseif((curva_piu_meno(x_VuExtra(ii),flag_extra,-1,velC_P).eq.0.d0).or.(curva_piu_meno(x_VuExtra(ii+1),flag_extra,-1,velC_P).eq.0.d0))then
                    iplog=2
                endif                                          
                
                position = ((rowIndex-1)*10)+ii+5
                position_d(ii+5)=position
                
                !!$omp parallel do
                DO ki=1,32
                    xtrasl(ki)=alfa*x_VuExtra(ki)+beta
                    A1(ki)=dmax1(0.d0,curva_piu_meno(xtrasl(ki),flag_extra,-1,velC_P))
                    B1(ki)=dmin1(estremo_m,curva_piu_meno(xtrasl(ki),flag_extra,1,velC_P))
                    alfa_j1(ki)=(B1(ki)-A1(ki))/2.d0
                    beta_j1(ki)=(B1(ki)+A1(ki))/2.d0                    
                    IF ((B1(ki)-A1(ki)).gt.10.d-14) THEN                                            
                        xinttrasl(ki)=(x_VuExtra(ki)+1.d0)*0.5d0
                        s(ki)=fi1(iplog,iqlog,xinttrasl(ki))
                        ds(ki)=dfi1(iplog,iqlog,xinttrasl(ki))
                    else            
                        xinttrasl(ki)=0.d0
                        s(ki)=0.d0
                        ds(ki)=0.d0
                    endif
                enddo
                !!$omp end parallel do

                !print *, rowIndex, ii, 'copy ep'
                ! istat = cudaMemcpy2DAsync(xtrasl_d(1,position),32,xtrasl(1),32,32,1,1,stream(position))
                ! istat = cudaMemcpy2DAsync(A1_d(1,position),32,A1(1),32,32,1,1,stream(position))
                ! istat = cudaMemcpy2DAsync(B1_d(1,position),32,B1(1),32,32,1,1,stream(position))
                ! istat = cudaMemcpy2DAsync(alfa_j1_d(1,position),32,alfa_j1(1),32,32,1,1,stream(position))
                ! istat = cudaMemcpy2DAsync(beta_j1_d(1,position),32,beta_j1(1),32,32,1,1,stream(position))
                ! istat = cudaMemcpy2DAsync(xinttrasl_d(1,position),32,xtrasl(1),32,32,1,1,stream(position))
                ! istat = cudaMemcpy2DAsync(s_d(1,position),32,s(1),32,32,1,1,stream(position))
                ! istat = cudaMemcpy2DAsync(ds_d(1,position),32,ds(1),32,32,1,1,stream(position))
                
                ! call performCalcEP<<<dimGrid,dimBlockCalculation,sharedMemDimension,stream(position)>>>(alfa_d(ii+5),position_d(ii+5),rowIndex_d)
                ! call sum_p2a<<<dimGrid,dimBlockCalculation,sharedMemDimension,stream(position)>>>(alfa_d(ii+5),position_d(ii+5), CalculationResults(ii+5,rowIndex))
                istat = cudaMemcpy2DAsync(xtrasl_d(1,position),32,xtrasl(1),32,32,1,1)
                istat = cudaMemcpy2DAsync(A1_d(1,position),32,A1(1),32,32,1,1)
                istat = cudaMemcpy2DAsync(B1_d(1,position),32,B1(1),32,32,1,1)
                istat = cudaMemcpy2DAsync(alfa_j1_d(1,position),32,alfa_j1(1),32,32,1,1)
                istat = cudaMemcpy2DAsync(beta_j1_d(1,position),32,beta_j1(1),32,32,1,1)
                istat = cudaMemcpy2DAsync(xinttrasl_d(1,position),32,xtrasl(1),32,32,1,1)
                istat = cudaMemcpy2DAsync(s_d(1,position),32,s(1),32,32,1,1)
                istat = cudaMemcpy2DAsync(ds_d(1,position),32,ds(1),32,32,1,1)
                
                call performCalcEP<<<dimGrid,dimBlockCalculation,sharedMemDimension>>>(alfa_d(ii+5),position_d(ii+5),rowIndex_d)
                call sum_p2a<<<dimGrid,dimBlockCalculation,sharedMemDimension>>>(alfa_d(ii+5),position_d(ii+5), CalculationResults(ii+5,rowIndex))
            endif
        enddo
        ! do ii=1,5
        !     if(xx(ii+1)-xx(ii).gt.1.d-14)then
        !         position = ((rowIndex-1)*10)+ii+5
        !         call performCalcEP<<<dimGrid,dimBlockCalculation,sharedMemDimension,stream(position)>>>(alfa_d(ii+5),position_d(ii+5),rowIndex_d)
        !         call sum_p2a<<<dimGrid,dimBlockCalculation,sharedMemDimension,stream(position)>>>(alfa_d(ii+5),position_d(ii+5), CalculationResults(ii+5,rowIndex_d))            
        !         ! ierrSync = cudaGetLastError()
        !         ! ierrAsync = cudaDeviceSynchronize()
        !         ! if (ierrSync /= cudaSuccess) write(*,*) 'Sync kernel error:', cudaGetErrorString(ierrSync)
        !         ! if (ierrAsync /= cudaSuccess) write(*,*) 'Async kernel error:', cudaGetErrorString(ierrAsync)
            
        !         ! app = CalculationResults(ii+5,rowIndex)
        !         ! print *, app, ii, rowIndex
        !     endif
        ! enddo
        !print *, 'performCalcEP'
        ! pause
  endif  
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