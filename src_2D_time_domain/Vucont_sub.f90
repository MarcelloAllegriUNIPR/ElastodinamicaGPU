Double precision function Vucont_sub(l_m_tilde,l_m,e_m_tilde,e_m,tempo1,tempo2,indice_i,indice_j)
!
!     INTEGRALE DELLE FUNZIONI DI FORMA li,lj SULL'ELEMENTO i
!     ELEMENTI SOVRAPPOSTI
!
  USE variable_2Dgeneral
  
  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABILI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  INTEGER(kind=4),INTENT(IN):: e_m_tilde, e_m, l_m_tilde, l_m, indice_i, indice_j
  REAL(kind=8),INTENT(IN):: tempo1, tempo2
  
  !Variabili locali
  INTEGER(kind=4) :: ind_gauss, Ngauss, AllocateStatus, iplog, iqlog, ki, kj, ip, iq, Ngaussi, ciclo, p, q
   
  REAL(kind=8) :: delta_x, xi1, xi2, xm1, xm2, xm, xn, alfa, alfa2, beta, Vucont_sub_S, Vucont_sub_P, int_1, int_2, & 
    p2, A1, B1, alfa_j1, beta_j1, p2a, s, ds, xtrasl, xinttrasl, r2_1, p3, p3a, cs1, serv, sm1, sm2, sn, sm, cs, cp
  
  REAL(kind=8),DIMENSION(:),ALLOCATABLE :: x, w, xint, wint, ti, vi, ww1
  
  REAL(kind=8),EXTERNAL :: fi1, dfi1, fiU
  
  
  INTEGER(kind=4):: l_serv, l_var, e_serv, e_var
  REAL(kind=8):: alpha, estremo_m_Vucont, estremo_m_tilde_Vucont, estremo_serv, estremo_var, th, tk, theta, theta2, &
                 coeff_1, coeff_2, coeff_3, coeff_4, coeff_5, a, b, sS1, sS2, sP1, sP2, aa, bb, a_estr_esterno, &
				 b_estr_esterno, dphi, cos_alpha, l_max
  INTEGER(kind=4), DIMENSION(2,2) :: delta_kronecker
  REAL(kind=8), DIMENSION(2) :: sigma, r_Vucont, vett_serv, vett_var
  REAL(kind=8), DIMENSION(2,10) :: ascisse_integrazione_S, ascisse_integrazione_SP, &
                                   ascisse_integrazione_S_1, ascisse_integrazione_S_2, &
                                   ascisse_integrazione_SP_1, ascisse_integrazione_SP_2								   
  INTEGER(kind=4), DIMENSION(2,10) :: ordinate_integrazione_S, ordinate_integrazione_SP, &
                                      ordinate_integrazione_S_1, ordinate_integrazione_S_2, &
									  ordinate_integrazione_SP_1, ordinate_integrazione_SP_2
                                         
  !!!!!!!!!!!!!!!!!!!!!!!!!! CORPO della SUBROUTINE !!!!!!!!!!!!!!!!!!
  cs=velC_S
  cp=velC_P
  
  estremo_m_Vucont=list_elements(e_m)%length
  estremo_m_tilde_Vucont=list_elements(e_m_tilde)%length
  
  if(e_m.lt.e_m_tilde)then
  !(e_m.lt.e_m_tilde.and.&
  !     (e_m.ne.1.and.e_m_tilde.ne.number_elements).or.&
	!   (e_m.eq.number_elements.and.e_m_tilde.eq.1))then
	 estremo_serv=estremo_m_Vucont
	 estremo_var=estremo_m_tilde_Vucont
	 l_serv=l_m
	 l_var=l_m_tilde
	 e_serv=e_m
     e_var=e_m_tilde
  else
	 estremo_serv=estremo_m_tilde_Vucont
	 estremo_var=estremo_m_Vucont
	 l_serv=l_m_tilde
	 l_var=l_m 
	 e_serv=e_m_tilde
     e_var=e_m
  endif

  
  theta=DACOS((list_nodes(e_serv+1)%coordinates(1)-list_nodes(e_serv)%coordinates(1))/estremo_serv)
  if(list_nodes(e_serv+1)%coordinates(2).lt.list_nodes(e_serv)%coordinates(2))then
     theta=-theta
  endif
 
  theta2=DACOS((list_nodes(e_var+1)%coordinates(1)-list_nodes(e_var)%coordinates(1))/estremo_var)
  if(list_nodes(e_var+1)%coordinates(2).lt.list_nodes(e_var)%coordinates(2))then
     theta2=-theta2
  endif
 
  vett_serv(1)=list_nodes(e_serv)%coordinates(1)-list_nodes(e_serv+1)%coordinates(1)
  vett_serv(2)=list_nodes(e_serv)%coordinates(2)-list_nodes(e_serv+1)%coordinates(2)
  vett_var(1)=list_nodes(e_var+1)%coordinates(1)-list_nodes(e_var)%coordinates(1)
  vett_var(2)=list_nodes(e_var+1)%coordinates(2)-list_nodes(e_var)%coordinates(2)
  alpha=DACOS((vett_serv(1)*vett_var(1)+vett_serv(2)*vett_var(2))/(estremo_serv*estremo_var)) 
  cos_alpha=(vett_serv(1)*vett_var(1)+vett_serv(2)*vett_var(2))/(estremo_serv*estremo_var)
  
  delta_kronecker(1,1)=1
  delta_kronecker(1,2)=0
  delta_kronecker(2,1)=0
  delta_kronecker(2,2)=1
  sigma(1)=DCOS(theta)
  sigma(2)=DSIN(theta)

	 ! write(*,*) 'alpha', alpha
     ! write(*,*) 'theta', theta
	 ! write(*,*) 'sigma1', sigma(1)
	 ! write(*,*) 'sigma2', sigma(2)

  coeff_1=((cp**2-cs**2)/(cp*cs))*(sigma(indice_i)*sigma(indice_j)-delta_kronecker(indice_i,indice_j)/2.d0)
  coeff_2=-((cp**2+cs**2)/(cp**2*cs**2))*delta_kronecker(indice_i,indice_j)/2.d0
  coeff_3=delta_kronecker(indice_i,indice_j)/2.d0
  coeff_4=sigma(indice_i)*sigma(indice_j)-delta_kronecker(indice_i,indice_j)/2.d0
  coeff_5=delta_kronecker(indice_i,indice_j)/2.d0
  
  Vucont_sub=0.d0
  delta_x=tempo1-tempo2
  IF (delta_x.le.0.d0) RETURN

  Vucont_sub_S=0.d0
  Vucont_sub_P=0.d0
  int_1=0.d0
  int_2=0.d0
   
  ascisse_integrazione_S=0.d0
  ordinate_integrazione_S=0
  ascisse_integrazione_SP=0.d0
  ordinate_integrazione_SP=0
  
  ascisse_integrazione_S_1=0.d0
  ascisse_integrazione_S_2=0.d0
  ascisse_integrazione_SP_1=0.d0
  ascisse_integrazione_SP_2=0.d0
  ordinate_integrazione_S_1=0
  ordinate_integrazione_S_2=0
  ordinate_integrazione_SP_1=0
  ordinate_integrazione_SP_2=0

  ind_gauss=6
  Ngauss=2**(ind_gauss-1) !128 nodi di Gauss

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

  ! Ngaussi=grado_q+1
  ! ALLOCATE(ti(Ngaussi),STAT=AllocateStatus)
  ! IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
  ! ALLOCATE(vi(Ngaussi),STAT=AllocateStatus)
  ! IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
  ! ti=gauss(Ngaussi)%nodiquad
  ! vi=gauss(Ngaussi)%pesiquad
  ! ALLOCATE(ww1(Ngaussi),STAT=AllocateStatus)
  ! IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
  
  Ngaussi=2**(grado_q)
  !Ngaussi=2**(6)
  ALLOCATE(ti(Ngaussi),STAT=AllocateStatus)
  IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
  ALLOCATE(vi(Ngaussi),STAT=AllocateStatus)
  IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
  ti=gauss(grado_q+1)%nodiquad
  vi=gauss(grado_q+1)%pesiquad
  !ti=gauss(6+1)%nodiquad
  !vi=gauss(6+1)%pesiquad
  ALLOCATE(ww1(Ngaussi),STAT=AllocateStatus)
  IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

  iplog=1
  iqlog=1
  
  !*************************************
  !             elementi allineati
  !*************************************
  
  if(dabs(cos_alpha+1).le.(10.d-6))then !se l'angolo alpha è circa pigreco ho elem. allineati
	 ! write(*,*) 'sono negli allineati'
	 ! write(*,*) 'elementi', e_m, e_m_tilde
	 ! write(*,*) 'alpha', alpha
	 ! write(*,*) 'cosbeta1', cos(theta)
	 ! write(*,*) 'sinbeta1', sin(theta)
	 ! pause
	 a=dmin1(cp*delta_x-estremo_serv,cs*delta_x)
	 b=dmax1(cp*delta_x-estremo_serv,cs*delta_x)
	 if(cp*delta_x-estremo_serv.le.0.d0)then
	     ascisse_integrazione_S(1,1)=0
	     ascisse_integrazione_S(2,1)=cs*delta_x
	     ordinate_integrazione_S(1,1)=2 !l'ordinata può essere 0, estremo_serv, 1 (retta_P) o 2(retta_S)  
	     ordinate_integrazione_S(2,1)=3
		  
	     ascisse_integrazione_SP(1,1)=0
	     ascisse_integrazione_SP(2,1)=cs*delta_x
	     ordinate_integrazione_SP(1,1)=1
	     ordinate_integrazione_SP(2,1)=2
		  
    	 ascisse_integrazione_SP(1,2)=cs*delta_x
	     ascisse_integrazione_SP(2,2)=cp*delta_x
	     ordinate_integrazione_SP(1,2)=1
	     ordinate_integrazione_SP(2,2)=3
		  
	 elseif(cs*delta_x-estremo_serv.gt.0.d0)then
         ascisse_integrazione_S(1,1)=0
	     ascisse_integrazione_S(2,1)=cs*delta_x-estremo_serv
	     ordinate_integrazione_S(1,1)=0   
	     ordinate_integrazione_S(2,1)=3
		  
	     ascisse_integrazione_S(1,2)=cs*delta_x-estremo_serv
	     ascisse_integrazione_S(2,2)=cs*delta_x
	     ordinate_integrazione_S(1,2)=2
	     ordinate_integrazione_S(2,2)=3
		  
	     ascisse_integrazione_SP(1,1)=cs*delta_x-estremo_serv
	     ascisse_integrazione_SP(2,1)=a
	     ordinate_integrazione_SP(1,1)=0
	     ordinate_integrazione_SP(2,1)=2
		  
	     ascisse_integrazione_SP(1,2)=a
	     ascisse_integrazione_SP(2,2)=cs*delta_x
	     ordinate_integrazione_SP(1,2)=1
	     ordinate_integrazione_SP(2,2)=2
		  
	     ascisse_integrazione_SP(1,3)=cs*delta_x
	     ascisse_integrazione_SP(2,3)=b
	     ordinate_integrazione_SP(1,3)=0
	     ordinate_integrazione_SP(2,3)=3
		  
         ascisse_integrazione_SP(1,4)=b
	     ascisse_integrazione_SP(2,4)=cp*delta_x
	     ordinate_integrazione_SP(1,4)=1
	     ordinate_integrazione_SP(2,4)=3
	   
     else
	     ascisse_integrazione_S(1,1)=0
	     ascisse_integrazione_S(2,1)=cs*delta_x
	     ordinate_integrazione_S(1,1)=2
	     ordinate_integrazione_S(2,1)=3
	  
	     ascisse_integrazione_SP(1,1)=0
	     ascisse_integrazione_SP(2,1)=a
	     ordinate_integrazione_SP(1,1)=0
	     ordinate_integrazione_SP(2,1)=2
		  
	     ascisse_integrazione_SP(1,2)=a
	     ascisse_integrazione_SP(2,2)=cs*delta_x
	     ordinate_integrazione_SP(1,2)=1
	     ordinate_integrazione_SP(2,2)=2
		  
	     ascisse_integrazione_SP(1,3)=cs*delta_x
	     ascisse_integrazione_SP(2,3)=b
	     ordinate_integrazione_SP(1,3)=0
	     ordinate_integrazione_SP(2,3)=3
		  
	     ascisse_integrazione_SP(1,4)=b
	     ascisse_integrazione_SP(2,4)=cp*delta_x
	     ordinate_integrazione_SP(1,4)=1
	     ordinate_integrazione_SP(2,4)=3
	   
     endif
	 
	 !!!!!!!!!!!!!!!!!!!!ciclo integrazione su insieme ES!!!!!!!!!!!!!!
	 DO ciclo=1,10 !10 è la dimensione massima degli array ascisse.../ordinate...	
  	
         if(ascisse_integrazione_S(1,ciclo).ge.estremo_var)then
	         Vucont_sub_S=Vucont_sub_S+0.d0
         elseif(ascisse_integrazione_S(1,ciclo).eq.ascisse_integrazione_S(2,ciclo))then
	         Vucont_sub_S=Vucont_sub_S+0.d0
         else
	         if(ascisse_integrazione_S(2,ciclo).ge.estremo_var)then
	             ascisse_integrazione_S(2,ciclo)=estremo_var
	         endif
	         
			 a_estr_esterno=ascisse_integrazione_S(1,ciclo)
			 b_estr_esterno=ascisse_integrazione_S(2,ciclo)
	         alfa=(b_estr_esterno-a_estr_esterno)/2.d0
	         beta=(b_estr_esterno+a_estr_esterno)/2.d0

             !secondo integrale   
             !------------------
	
             p2 = 0.d0
			 
			 if(ordinate_integrazione_S(1,ciclo).ne.0)then
			     iplog=2
				 iqlog=1
			 endif
	 
    	     DO ki=1,Ngauss 
	
                 xtrasl=alfa*x(ki)+beta
			 
		    	 A1=retta_S_P(xtrasl,ordinate_integrazione_S(1,ciclo))
           	     B1=retta_S_P(xtrasl,ordinate_integrazione_S(2,ciclo))
			 
	             alfa_j1=(B1-A1)/2.d0
	             beta_j1=(B1+A1)/2.d0
                  
				 p2a = 0.d0
	             IF ((B1-A1).gt.10.d-15) THEN
	
                     DO kj=1,Ngauss

	                     xinttrasl=(xint(kj)+1.d0)*0.5d0
	                     s=fi1(iplog,iqlog,xinttrasl)
	                     ds=dfi1(iplog,iqlog,xinttrasl)
	                     serv=alfa_j1*(2.d0*s-1.d0)+beta_j1
                         r2_1=(estremo_serv+xtrasl-serv)**2

	                     p2a = p2a + wint(kj)*ds*fiU(l_serv,serv,estremo_serv,grado_q)*coeff_1*delta_x*(cp*dsqrt((cs**2)*(delta_x**2)-r2_1)+cs*dsqrt((cp**2)*(delta_x**2)-r2_1))**(-1)
	    
	                     IF (delta_kronecker(indice_i,indice_j).eq.1) THEN
	                         ! write(*,*) 'indice croneker', indice_i, indice
	                         p2a= p2a+wint(kj)*ds*fiU(l_serv,serv,estremo_serv,grado_q)*coeff_3*((1.d0/(cp**2))*dlog(cp*delta_x+dsqrt(dabs(-r2_1+cp**2*delta_x**2)))+(1.d0/(cs**2))*dlog(cs*delta_x+dsqrt(dabs(-r2_1+cs**2*delta_x**2))))
	                     ENDIF
		  
                     END DO

        	     ENDIF

	             p2 = p2 +p2a*alfa_j1*w(ki)*fiU(l_var,xtrasl,estremo_var,grado_q)

             END DO
             Vucont_sub_S=Vucont_sub_S+p2*alfa
			 ! write(*,*) 'elementi', e_m, e_m_tilde
             ! write(*,*) 'p2', p2 			 
  
               !primo integrale     
               !-----------------

	         p3 = 0.d0

			 if(a_estr_esterno.eq.0.d0)then
			     p=2
				 q=1
			 endif			 
    
             IF (delta_kronecker(indice_i,indice_j).eq.1) THEN
            ! write(*,*) 'indice croneker', indice_i, indice_j
                 DO ki=1,Ngauss
	    
	                 p3a = 0.d0
	 
	                 xtrasl=2.d0*alfa*fi1(p,q,(x(ki)+1)/2.d0)+a_estr_esterno
	                 dphi=dfi1(p,q,(x(ki)+1)/2.d0)
					 
				     A1=retta_S_P(xtrasl,ordinate_integrazione_S(1,ciclo))
           	         B1=retta_S_P(xtrasl,ordinate_integrazione_S(2,ciclo))

	                 alfa_j1=(B1-A1)/2.d0
	                 beta_j1=(B1+A1)/2.d0
	  
	                 if ((B1-A1).gt.10.d-15) then

		                 cs1=(estremo_serv+xtrasl-beta_j1)/alfa_j1

	                     if (dabs(cs1).lt.1.25d0) then

		                     call mopelog(Ngaussi,ti,cs1,0.d0,ww1,1)
    	 
                             DO kj=1,Ngaussi
			                     p3a=p3a+vi(kj)*ww1(kj)*fiU(l_serv,alfa_j1*ti(kj)+beta_j1,estremo_serv,grado_q)
                             END DO

	                         DO kj=1,Ngauss
	                             p3a=p3a+w(kj)*fiU(l_serv,alfa_j1*x(kj)+beta_j1,estremo_serv,grado_q)*dlog(alfa_j1)
                             END DO

	                     else
    	                     !Write(*,*) 'non ho usato la mopelog, allineati S'
				             !pause
	                         DO kj=1,Ngauss
	                             p3a=p3a+w(kj)*fiU(l_serv,alfa_j1*x(kj)+beta_j1,estremo_serv,grado_q)*(dlog(dabs(cs1-x(kj)))+dlog(alfa_j1))
                             END DO
                
	                     endif

	                 endif

                     p3=p3+p3a*alfa_j1*w(ki)*fiU(l_var,xtrasl,estremo_var,grado_q)*dphi
                 END DO
    
             ENDIF  
             Vucont_sub_S=Vucont_sub_S+coeff_2*p3*alfa
		 
	     ENDIF
     END DO  
  
     !!!!!!!!!!!!!!!!!!!!ciclo integrazione su insieme ESP!!!!!!!!!!!!!!
     DO ciclo=1,10 !10 è la dimensione massima degli array ascisse.../ordinate...	
  	
         if(ascisse_integrazione_SP(1,ciclo).ge.estremo_var)then
	         Vucont_sub_P=Vucont_sub_P+0.d0
         elseif(ascisse_integrazione_SP(1,ciclo).eq.ascisse_integrazione_SP(2,ciclo))then
	         Vucont_sub_P=Vucont_sub_P+0.d0
         else
	         if(ascisse_integrazione_SP(2,ciclo).ge.estremo_var)then
	             ascisse_integrazione_SP(2,ciclo)=estremo_var
	         endif
	
	         alfa=(ascisse_integrazione_SP(2,ciclo)-ascisse_integrazione_SP(1,ciclo))/2.d0
	         beta=(ascisse_integrazione_SP(2,ciclo)+ascisse_integrazione_SP(1,ciclo))/2.d0

             !secondo integrale   
             !------------------
	
             p2 = 0.d0
	         
			 if(ordinate_integrazione_SP(1,ciclo).ne.0)then
			     iplog=2
				 iqlog=1
			 endif
	 
	         DO ki=1,Ngauss 
	
                 xtrasl=alfa*x(ki)+beta
			 
			     A1=retta_S_P(xtrasl,ordinate_integrazione_SP(1,ciclo))
           	     B1=retta_S_P(xtrasl,ordinate_integrazione_SP(2,ciclo))
			 
	             alfa_j1=(B1-A1)/2.d0
	             beta_j1=(B1+A1)/2.d0

				 p2a = 0.d0
	             IF ((B1-A1).gt.10.d-15) THEN
	
	                
                     DO kj=1,Ngauss

	                     xinttrasl=(xint(kj)+1.d0)*0.5d0
	                     s=fi1(iplog,iqlog,xinttrasl)
	                     ds=dfi1(iplog,iqlog,xinttrasl)
	                     serv=alfa_j1*(2.d0*s-1.d0)+beta_j1
                         r2_1=(estremo_serv+xtrasl-serv)**2

	                     p2a = p2a + wint(kj)*ds*fiU(l_serv,serv,estremo_serv,grado_q)*coeff_4*(((cp**2)*r2_1)**(-1))*cp*delta_x*dsqrt(dabs(-r2_1+cp**2*delta_x**2))
	    
	                     IF (delta_kronecker(indice_i,indice_j).eq.1) THEN
	    
	                         p2a= p2a+wint(kj)*ds*fiU(l_serv,serv,estremo_serv,grado_q)*coeff_5*(1.d0/(cp**2))*dlog(cp*delta_x+dsqrt(dabs(-r2_1+cp**2*delta_x**2)))
	                     ENDIF
		  
                     END DO
 
            	 ENDIF

	             p2 = p2 +p2a*alfa_j1*w(ki)*fiU(l_var,xtrasl,estremo_var,grado_q)

             END DO
             Vucont_sub_P=Vucont_sub_P+p2*alfa
  
             ! primo integrale     
             !-----------------

	         p3 = 0.d0
    
             IF (delta_kronecker(indice_i,indice_j).eq.1) THEN
    
                 DO ki=1,Ngauss
	    
	                 p3a = 0.d0
	 
	                 xtrasl=alfa*x(ki)+beta
	             
				     A1=retta_S_P(xtrasl,ordinate_integrazione_SP(1,ciclo))
           	         B1=retta_S_P(xtrasl,ordinate_integrazione_SP(2,ciclo))
			 
	                 alfa_j1=(B1-A1)/2.d0
	                 beta_j1=(B1+A1)/2.d0
	  
	                 if ((B1-A1).gt.10.d-15) then

		                 cs1=(estremo_serv+xtrasl-beta_j1)/alfa_j1

	                     if (dabs(cs1).lt.1.25d0) then

		                     call mopelog(Ngaussi,ti,cs1,0.d0,ww1,1)
    	 
                             DO kj=1,Ngaussi
			                     p3a=p3a+vi(kj)*ww1(kj)*fiU(l_serv,alfa_j1*ti(kj)+beta_j1,estremo_serv,grado_q)
                             END DO

	                         DO kj=1,Ngauss
	                             p3a=p3a+w(kj)*fiU(l_serv,alfa_j1*x(kj)+beta_j1,estremo_serv,grado_q)*dlog(alfa_j1)
                             END DO

	                     else
    	                     !Write(*,*) 'non ho usato la mopelog, allineati P'
				             !pause
	                         DO kj=1,Ngauss
	                             p3a=p3a+w(kj)*fiU(l_serv,alfa_j1*x(kj)+beta_j1,estremo_serv,grado_q)*(dlog(dabs(cs1-x(kj)))+dlog(alfa_j1))
                             END DO
                
	                     endif

	                 endif

                     p3=p3+p3a*alfa_j1*w(ki)*fiU(l_var,xtrasl,estremo_var,grado_q)
                 END DO
    
             ENDIF  
             Vucont_sub_P=Vucont_sub_P-(1.d0/(cp**2))*coeff_5*p3*alfa
		 
	     ENDIF
     END DO
	 ! write(*,*) 'integrale in S', Vucont_sub_S
	 ! write(*,*) 'integrale in P', Vucont_sub_P
	 ! pause
	 
  ! !*************************************
  ! !             elementi disallineati 
  ! !************************************* 
  
  else
     !write(*,*) 'elementi', e_m, e_m_tilde
     !write(*,*) 'alpha', alpha
     !write(*,*) 'sono nei disallineati'
	 !pause
     if((-4*(estremo_serv*SIN(alpha))**2+4*(cs*delta_x)**2).ge.0.d0)then
	     sS1=estremo_serv*COS(alpha)-dsqrt(-(estremo_serv*SIN(alpha))**2+(cs*delta_x)**2)
         sS2=estremo_serv*COS(alpha)+dsqrt(-(estremo_serv*SIN(alpha))**2+(cs*delta_x)**2)		 
	 else
	     sS1=0.d0
		 sS2=0.d0
	 endif
	 
	 if((-4*(estremo_serv*SIN(alpha))**2+4*(cp*delta_x)**2).ge.0.d0)then
	     sP1=estremo_serv*COS(alpha)-dsqrt(-(estremo_serv*SIN(alpha))**2+(cp*delta_x)**2)
         sP2=estremo_serv*COS(alpha)+dsqrt(-(estremo_serv*SIN(alpha))**2+(cp*delta_x)**2)		 
	 else
	     sP1=0.d0
		 sP2=0.d0
	 endif
	 
	 a=dmin1(cp*delta_x,cs*delta_x/dabs(SIN(alpha)))
	 b=dmax1(cp*delta_x,cs*delta_x/dabs(SIN(alpha)))
	 
	 if((alpha.gt.0.d0.AND.alpha.lt.(8.d0*datan(1.d0)/4.d0)).OR.(alpha.gt.(3.d0*8.d0*datan(1.d0)/4.d0).AND.alpha.lt.(8.d0*datan(1.d0))))then
	     !write(*,*) 'sono nei disallineati caso sbagliato'
		 !write(*,*) 'elementi', e_m, e_m_tilde
		 !write(*,*) 'alpha', alpha
		 ! write(*,*) 'cosbeta1', cos(theta)
		 ! write(*,*) 'sinbeta1', sin(theta)
		 ! write(*,*) 'cosbeta2', cos(theta2)
		 ! write(*,*) 'sinbeta2', sin(theta2)
		 !pause
		 
		 l_max=dmax1(estremo_serv,estremo_var)
		 if(cs*delta_x.ge.l_max.and.(estremo_serv**2+estremo_var**2&
		             -2.d0*cos(alpha)*estremo_serv*estremo_var.le.(cs*delta_x)**2))then
		 ascisse_integrazione_S_1(1,1)=0.d0
		 ascisse_integrazione_S_1(2,1)=estremo_var
		 ordinate_integrazione_S_1(1,1)=0
		 ordinate_integrazione_S_1(2,1)=5 
		
		 else
		 ascisse_integrazione_S_1(1,1)=0.d0
		 ascisse_integrazione_S_1(2,1)=cs*delta_x
		 ordinate_integrazione_S_1(1,1)=2
		 ordinate_integrazione_S_1(2,1)=5
		 
		 ascisse_integrazione_S_1(1,2)=cs*delta_x
		 ascisse_integrazione_S_1(2,2)=(cs*delta_x)/dabs(SIN(alpha))
		 ordinate_integrazione_S_1(1,2)=2
		 ordinate_integrazione_S_1(2,2)=3
		 
		 if(estremo_serv-(cs*delta_x)/dabs(SIN(alpha)).lt.0.d0)then
		     if(estremo_serv-cs*delta_x*COS(alpha)/dabs(SIN(alpha)).gt.0.d0)then
		         ascisse_integrazione_S_2(1,1)=dmax1(0.d0,sS1)
			     ascisse_integrazione_S_2(2,1)=sS2
			     ordinate_integrazione_S_2(1,1)=2
			     ordinate_integrazione_S_2(2,1)=0
		     else	 
			     ascisse_integrazione_S_2(1,1)=dmax1(0.d0,sS1)
			     ascisse_integrazione_S_2(2,1)=sS2
			     ordinate_integrazione_S_2(1,1)=2
			     ordinate_integrazione_S_2(2,1)=0
			 
			     ascisse_integrazione_S_2(1,2)=sS2
			     ascisse_integrazione_S_2(2,2)=(cs*delta_x)/dabs(SIN(alpha))
			     ordinate_integrazione_S_2(1,2)=2
			     ordinate_integrazione_S_2(2,2)=3
		     endif
			 
	     endif
		 
		 ascisse_integrazione_SP_1(1,1)=0.d0
		 ascisse_integrazione_SP_1(2,1)=cs*delta_x/dabs(SIN(alpha))
		 ordinate_integrazione_SP_1(1,1)=1
		 ordinate_integrazione_SP_1(2,1)=2
		 
		 ascisse_integrazione_SP_1(1,2)=cs*delta_x
		 ascisse_integrazione_SP_1(2,2)=a
		 ordinate_integrazione_SP_1(1,2)=3
		 ordinate_integrazione_SP_1(2,2)=5
		 
		 ascisse_integrazione_SP_1(1,3)=a
		 ascisse_integrazione_SP_1(2,3)=cp*delta_x
		 ordinate_integrazione_SP_1(1,3)=1
		 ordinate_integrazione_SP_1(2,3)=5
		 
		 ascisse_integrazione_SP_1(1,4)=cp*delta_x
		 ascisse_integrazione_SP_1(2,4)=b
		 ordinate_integrazione_SP_1(1,4)=3
		 ordinate_integrazione_SP_1(2,4)=4

		 ascisse_integrazione_SP_1(1,5)=b
		 ascisse_integrazione_SP_1(2,5)=cp*delta_x/dabs(SIN(alpha))
		 ordinate_integrazione_SP_1(1,5)=1
		 ordinate_integrazione_SP_1(2,5)=4
		
		 
		 if(estremo_serv-(cp*delta_x)/dabs(SIN(alpha)).lt.0.d0)then
		     aa=dmin1(estremo_serv-cs*delta_x/dabs(SIN(alpha)),estremo_serv-cp*delta_x*COS(alpha)/dabs(SIN(alpha)))
			 bb=dmax1(estremo_serv-cs*delta_x/dabs(SIN(alpha)),estremo_serv-cp*delta_x*COS(alpha)/dabs(SIN(alpha)))
			 
		     if(aa.gt.0.d0)then
			 
     			 ascisse_integrazione_SP_2(1,1)=dmax1(0.d0,sP1)
			     ascisse_integrazione_SP_2(2,1)=sP2
			     ordinate_integrazione_SP_2(1,1)=1
			     ordinate_integrazione_SP_2(2,1)=0
		     
			 elseif((aa.le.0.d0).AND.(bb.ge.0.d0))then 
			     if(bb==estremo_serv-cs*delta_x/dabs(SIN(alpha)))then
				     
					 ascisse_integrazione_SP_2(1,1)=dmax1(0.d0,sP1)
			         ascisse_integrazione_SP_2(2,1)=sP2
			         ordinate_integrazione_SP_2(1,1)=1
			         ordinate_integrazione_SP_2(2,1)=0
					 
					 ascisse_integrazione_SP_2(1,2)=sP2
			         ascisse_integrazione_SP_2(2,2)=cp*delta_x/dabs(SIN(alpha))
			         ordinate_integrazione_SP_2(1,2)=1
			         ordinate_integrazione_SP_2(2,2)=4
			     else
				 
				     ascisse_integrazione_SP_2(1,1)=dmax1(0.d0,sP1)
			         ascisse_integrazione_SP_2(2,1)=dmax1(0.d0,sS1)
			         ordinate_integrazione_SP_2(1,1)=1
			         ordinate_integrazione_SP_2(2,1)=0
					 
					 ascisse_integrazione_SP_2(1,2)=dmax1(0.d0,sS1)
			         ascisse_integrazione_SP_2(2,2)=sS2
			         ordinate_integrazione_SP_2(1,2)=1
			         ordinate_integrazione_SP_2(2,2)=2
					 
					 ascisse_integrazione_SP_2(1,3)=sS2
			         ascisse_integrazione_SP_2(2,3)=sP2
			         ordinate_integrazione_SP_2(1,3)=1
			         ordinate_integrazione_SP_2(2,3)=0
	                 		     
				 endif
			 elseif((bb.lt.0.d0).AND.(estremo_serv-cs*delta_x*COS(alpha)/dabs(SIN(alpha)).gt.0.d0))then
                 
				 ascisse_integrazione_SP_2(1,1)=dmax1(0.d0,sP1)
			     ascisse_integrazione_SP_2(2,1)=dmax1(0.d0,sS1)
			     ordinate_integrazione_SP_2(1,1)=1
		         ordinate_integrazione_SP_2(2,1)=0
					 
	    		 ascisse_integrazione_SP_2(1,2)=dmax1(0.d0,sS1)
		         ascisse_integrazione_SP_2(2,2)=sS2
			     ordinate_integrazione_SP_2(1,2)=1
		         ordinate_integrazione_SP_2(2,2)=2
					 
				 ascisse_integrazione_SP_2(1,3)=sS2
	             ascisse_integrazione_SP_2(2,3)=sP2
		         ordinate_integrazione_SP_2(1,3)=1
		         ordinate_integrazione_SP_2(2,3)=0
				 
				 ascisse_integrazione_SP_2(1,4)=sP2
		         ascisse_integrazione_SP_2(2,4)=cp*delta_x/dabs(SIN(alpha))
			     ordinate_integrazione_SP_2(1,4)=1
		         ordinate_integrazione_SP_2(2,4)=4
				 
             else
				 ascisse_integrazione_SP_2(1,1)=dmax1(0.d0,sP1)
			     ascisse_integrazione_SP_2(2,1)=dmax1(0.d0,sS1)
			     ordinate_integrazione_SP_2(1,1)=1
		         ordinate_integrazione_SP_2(2,1)=0
					 
	    		 ascisse_integrazione_SP_2(1,2)=dmax1(0.d0,sS1)
		         ascisse_integrazione_SP_2(2,2)=cs*delta_x/dabs(SIN(alpha))
			     ordinate_integrazione_SP_2(1,2)=1
		         ordinate_integrazione_SP_2(2,2)=2
				 
				 ascisse_integrazione_SP_2(1,3)=sS2
			     ascisse_integrazione_SP_2(2,3)=dmin1(sP2,cs*delta_x/dabs(SIN(alpha)))
			     ordinate_integrazione_SP_2(1,3)=3
		         ordinate_integrazione_SP_2(2,3)=0
					 
	    		 ascisse_integrazione_SP_2(1,4)=dmin1(sP2,cs*delta_x/dabs(SIN(alpha)))
		         ascisse_integrazione_SP_2(2,4)=sP2
			     ordinate_integrazione_SP_2(1,4)=1
		         ordinate_integrazione_SP_2(2,4)=0
				 
				 ascisse_integrazione_SP_2(1,5)=sP2
			     ascisse_integrazione_SP_2(2,5)=dmax1(sP2,cs*delta_x/dabs(SIN(alpha)))
			     ordinate_integrazione_SP_2(1,5)=3
		         ordinate_integrazione_SP_2(2,5)=4
					 
	    		 ascisse_integrazione_SP_2(1,6)=dmax1(sP2,cs*delta_x/dabs(SIN(alpha)))
		         ascisse_integrazione_SP_2(2,6)=cp*delta_x/dabs(SIN(alpha))
			     ordinate_integrazione_SP_2(1,6)=1
		         ordinate_integrazione_SP_2(2,6)=4                 				 
				 
			 endif	 
		 endif
		 endif
		 !endif

		 
	 elseif(((alpha.ge.(8.d0*datan(1.d0)/4.d0)).AND.(alpha.lt.(8.d0*datan(1.d0)/2.d0))).OR.((alpha.gt.(8.d0*datan(1.d0)/2.d0)).AND.(alpha.le.(3.d0*8.d0*datan(1.d0)/4.d0))))then
	    !write(*,*) 'sono nei disallineati caso giusto'
		!write(*,*) 'elementi', e_m, e_m_tilde
		!write(*,*) 'alpha', alpha
		!write(*,*) 'cosbeta1', cos(theta)
		!write(*,*) 'sinbeta1', sin(theta)
		!pause
		 !write(*,*) 'angolo', alpha
		 
		 a=dmin1(sP2,cs*delta_x)
	     b=dmax1(sP2,cs*delta_x)
	     if(sP2.le.0.d0)then
	         ascisse_integrazione_S_1(1,1)=0.d0
	         ascisse_integrazione_S_1(2,1)=cs*delta_x
	         ordinate_integrazione_S_1(1,1)=2 !l'ordinata può essere 0, estremo_serv, 1 (retta_P) o 2(retta_S)  
	         ordinate_integrazione_S_1(2,1)=5
		  
	         ascisse_integrazione_SP_1(1,1)=0.d0
	         ascisse_integrazione_SP_1(2,1)=cs*delta_x
	         ordinate_integrazione_SP_1(1,1)=1
	         ordinate_integrazione_SP_1(2,1)=2
		  
    	     ascisse_integrazione_SP_1(1,2)=cs*delta_x
	         ascisse_integrazione_SP_1(2,2)=cp*delta_x
	         ordinate_integrazione_SP_1(1,2)=1
	         ordinate_integrazione_SP_1(2,2)=5
		  
	     elseif(sS2.gt.0.d0)then
             ascisse_integrazione_S_1(1,1)=0.d0
	         ascisse_integrazione_S_1(2,1)=sS2
	         ordinate_integrazione_S_1(1,1)=0   
	         ordinate_integrazione_S_1(2,1)=5
		  
	         ascisse_integrazione_S_1(1,2)=sS2
	         ascisse_integrazione_S_1(2,2)=cs*delta_x
	         ordinate_integrazione_S_1(1,2)=2
	         ordinate_integrazione_S_1(2,2)=5
		  
	         ascisse_integrazione_SP_1(1,1)=sS2
	         ascisse_integrazione_SP_1(2,1)=a
	         ordinate_integrazione_SP_1(1,1)=0
	         ordinate_integrazione_SP_1(2,1)=2
		  
	         ascisse_integrazione_SP_1(1,2)=a
	         ascisse_integrazione_SP_1(2,2)=cs*delta_x
	         ordinate_integrazione_SP_1(1,2)=1
	         ordinate_integrazione_SP_1(2,2)=2
		  
	         ascisse_integrazione_SP_1(1,3)=cs*delta_x
	         ascisse_integrazione_SP_1(2,3)=b
	         ordinate_integrazione_SP_1(1,3)=0
	         ordinate_integrazione_SP_1(2,3)=5
		  
             ascisse_integrazione_SP_1(1,4)=b
	         ascisse_integrazione_SP_1(2,4)=cp*delta_x
	         ordinate_integrazione_SP_1(1,4)=1
	         ordinate_integrazione_SP_1(2,4)=5
	   
         else
	         ascisse_integrazione_S_1(1,1)=0.d0
	         ascisse_integrazione_S_1(2,1)=cs*delta_x
	         ordinate_integrazione_S_1(1,1)=2
	         ordinate_integrazione_S_1(2,1)=5
	  
	         ascisse_integrazione_SP_1(1,1)=0.d0
	         ascisse_integrazione_SP_1(2,1)=a
	         ordinate_integrazione_SP_1(1,1)=0
	         ordinate_integrazione_SP_1(2,1)=2
		  
	         ascisse_integrazione_SP_1(1,2)=a
	         ascisse_integrazione_SP_1(2,2)=cs*delta_x
	         ordinate_integrazione_SP_1(1,2)=1
	         ordinate_integrazione_SP_1(2,2)=2
		  
	         ascisse_integrazione_SP_1(1,3)=cs*delta_x
	         ascisse_integrazione_SP_1(2,3)=b
	         ordinate_integrazione_SP_1(1,3)=0
	         ordinate_integrazione_SP_1(2,3)=5
		  
	         ascisse_integrazione_SP_1(1,4)=b
	         ascisse_integrazione_SP_1(2,4)=cp*delta_x
	         ordinate_integrazione_SP_1(1,4)=1
	         ordinate_integrazione_SP_1(2,4)=5
	   
         endif
		 ! ascisse_integrazione_S_1(1,1)=0.d0
		 ! ascisse_integrazione_S_1(2,1)=cs*delta_x
		 ! ordinate_integrazione_S_1(1,1)=2
		 ! ordinate_integrazione_S_1(2,1)=5
		 
		 ! if(estremo_serv-cs*delta_x.lt.0.d0)then
		     ! ascisse_integrazione_S_2(1,1)=0.d0
			 ! ascisse_integrazione_S_2(2,1)=sS2
			 ! ordinate_integrazione_S_2(1,1)=2
			 ! ordinate_integrazione_S_2(2,1)=0
		 ! endif
		 
		 ! ascisse_integrazione_SP_1(1,1)=0.d0
		 ! ascisse_integrazione_SP_1(2,1)=cs*delta_x
		 ! ordinate_integrazione_SP_1(1,1)=1
		 ! ordinate_integrazione_SP_1(2,1)=2
		 
		 ! ascisse_integrazione_SP_1(1,2)=cs*delta_x
		 ! ascisse_integrazione_SP_1(2,2)=cp*delta_x
		 ! ordinate_integrazione_SP_1(1,2)=1
		 ! ordinate_integrazione_SP_1(2,2)=5
		 
		 ! if(estremo_serv-cp*delta_x.lt.0.d0)then
             ! if((0.d0).lt.estremo_serv-cs*delta_x)then
				 ! ascisse_integrazione_SP_2(1,1)=0.d0
		         ! ascisse_integrazione_SP_2(2,1)=sP2
		         ! ordinate_integrazione_SP_2(1,1)=1
		         ! ordinate_integrazione_SP_2(2,1)=0
		     ! else
				 ! ascisse_integrazione_SP_2(1,1)=0.d0
		         ! ascisse_integrazione_SP_2(2,1)=sS2
		         ! ordinate_integrazione_SP_2(1,1)=1
		         ! ordinate_integrazione_SP_2(2,1)=2
				 
				 ! ascisse_integrazione_SP_2(1,2)=sS2
		         ! ascisse_integrazione_SP_2(2,2)=sP2
		         ! ordinate_integrazione_SP_2(1,2)=1
		         ! ordinate_integrazione_SP_2(2,2)=0	
             ! endif				 
		 ! endif
		 
		 
	 endif
	 
	 !!!!!!!!!!!!!!!!!!!!!!ciclo integrazione su insieme ES!!!!!!!!!!!!!!
     DO ciclo=1,10 !10 è la dimensione massima degli array ascisse.../ordinate...	
  	     if(ascisse_integrazione_S_1(1,ciclo).ge.estremo_var)then
	         int_1=int_1+0.d0
         elseif(ascisse_integrazione_S_1(1,ciclo) == ascisse_integrazione_S_1(2,ciclo))then
	         int_1=int_1+0.d0
         else
	         if(ascisse_integrazione_S_1(2,ciclo).ge.estremo_var)then
	             ascisse_integrazione_S_1(2,ciclo)=estremo_var
	         endif
	         
			 a_estr_esterno=ascisse_integrazione_S_1(1,ciclo)
			 b_estr_esterno=ascisse_integrazione_S_1(2,ciclo)
	         alfa=(b_estr_esterno-a_estr_esterno)/2.d0
	         beta=(b_estr_esterno+a_estr_esterno)/2.d0

             !secondo integrale   
             !------------------
	
             p2 = 0.d0
	         
			 if(ordinate_integrazione_S_1(1,ciclo).ne.0)then
			     iplog=3
				 iqlog=1
			 endif
             if(ordinate_integrazione_S_1(2,ciclo).ne.5)then
			     iplog=1
				 iqlog=3
			 endif
	         
			 DO ki=1,Ngauss 
	
                 xtrasl=alfa*x(ki)+beta
			 
			     A1=arco_S_P(xtrasl,ordinate_integrazione_S_1(1,ciclo))
           	     B1=arco_S_P(xtrasl,ordinate_integrazione_S_1(2,ciclo))
			 
	             alfa_j1=(B1-A1)/2.d0
	             beta_j1=(B1+A1)/2.d0

				 p2a = 0.d0
				 
	             IF ((B1-A1).gt.10.d-15) THEN

                     DO kj=1,Ngauss

	                     xinttrasl=(xint(kj)+1.d0)*0.5d0
	                     s=fi1(iplog,iqlog,xinttrasl)
	                     ds=dfi1(iplog,iqlog,xinttrasl)
	                     serv=alfa_j1*(2.d0*s-1.d0)+beta_j1
                         r2_1=(estremo_serv-serv)**2+xtrasl**2-2*(estremo_serv-serv)*xtrasl*COS(alpha)
						 r_Vucont(1)=(estremo_serv-serv)*COS(theta)+xtrasl*COS(theta2)
						 r_Vucont(2)=(estremo_serv-serv)*SIN(theta)+xtrasl*SIN(theta2)
		   
	                     p2a = p2a + wint(kj)*ds*fiU(l_serv,serv,estremo_serv,grado_q)*((cp**2-cs**2)/(cp*cs))*delta_x*(cp*dsqrt((cs**2)*(delta_x**2)-r2_1)+cs*dsqrt((cp**2)*(delta_x**2)-r2_1))**(-1)*(r_Vucont(indice_i)*r_Vucont(indice_j)/r2_1-delta_kronecker(indice_i,indice_j)/2.d0)
	    
	                     IF (delta_kronecker(indice_i,indice_j).eq.1.d0) THEN
	    
	                         p2a= p2a+wint(kj)*ds*fiU(l_serv,serv,estremo_serv,grado_q)*coeff_3*((1.d0/(cp**2))*dlog(cp*delta_x+dsqrt(dabs(-r2_1+cp**2*delta_x**2)))+(1.d0/(cs**2))*dlog(cs*delta_x+dsqrt(dabs(-r2_1+cs**2*delta_x**2))))
	                     ENDIF
		  
                     END DO

        	     ENDIF

	             p2 = p2 +p2a*alfa_j1*w(ki)*fiU(l_var,xtrasl,estremo_var,grado_q)

             END DO
             int_1=int_1+p2*alfa
  
             !primo integrale     
             !-----------------

	         p3 = 0.d0
    
             IF (delta_kronecker(indice_i,indice_j).eq.1.d0) THEN
                 if(a_estr_esterno.eq.0.d0)then
			         p=3
				     q=1
			     endif			
				 
                 DO ki=1,Ngauss
	     
	                 p3a = 0.d0
	 
	                 xtrasl=2.d0*alfa*fi1(p,q,(x(ki)+1)/2.d0)+a_estr_esterno
	                 dphi=dfi1(p,q,(x(ki)+1)/2.d0)	             
				 
				     A1=arco_S_P(xtrasl,ordinate_integrazione_S_1(1,ciclo))
           	         B1=arco_S_P(xtrasl,ordinate_integrazione_S_1(2,ciclo))

	                 alfa_j1=(B1-A1)/2.d0
	                 beta_j1=(B1+A1)/2.d0
	  
	                 if ((B1-A1).gt.10.d-15) then
					     
					     DO kj=1,Ngauss
						     
							 iplog=2
							 iqlog=1
					         xinttrasl=(xint(kj)+1.d0)*0.5d0
	                         s=fi1(iplog,iqlog,xinttrasl)
	                         ds=dfi1(iplog,iqlog,xinttrasl)
	                         serv=alfa_j1*(2.d0*s-1.d0)+beta_j1
                             r2_1=(estremo_serv-serv)**2+xtrasl**2-2*(estremo_serv-serv)*xtrasl*COS(alpha)
						     
							 p3a=p3a+wint(kj)*ds*fiU(l_serv,serv,estremo_serv,grado_q)*dlog(dsqrt(r2_1))
							 
						 END DO
						 
	                 endif

                     p3=p3+p3a*alfa_j1*w(ki)*fiU(l_var,xtrasl,estremo_var,grado_q)*dphi
                 END DO
    
             ENDIF  
             int_1=int_1+coeff_2*p3*alfa
		 
	     ENDIF
     END DO  

     DO ciclo=1,10 !10 è la dimensione massima degli array ascisse.../ordinate...	
  	     if(ascisse_integrazione_S_2(1,ciclo).ge.estremo_var)then
	         int_2=int_2+0.d0
         elseif(ascisse_integrazione_S_2(1,ciclo) == ascisse_integrazione_S_2(2,ciclo))then
	         int_2=int_2+0.d0
         else
	         if(ascisse_integrazione_S_2(2,ciclo).ge.estremo_var)then
	             ascisse_integrazione_S_2(2,ciclo)=estremo_var
	         endif
	 
			 a_estr_esterno=ascisse_integrazione_S_2(1,ciclo)
			 b_estr_esterno=ascisse_integrazione_S_2(2,ciclo)
	         alfa=(b_estr_esterno-a_estr_esterno)/2.d0
	         beta=(b_estr_esterno+a_estr_esterno)/2.d0

             !secondo integrale   
             !------------------
	
             p2 = 0.d0
	         
			 if(ordinate_integrazione_S_2(1,ciclo).ne.0)then
			     iplog=3
				 iqlog=1
			 endif
             if(ordinate_integrazione_S_2(2,ciclo).ne.5)then
			     iplog=1
				 iqlog=3
			 endif
			 
	         DO ki=1,Ngauss 
	
                 xtrasl=alfa*x(ki)+beta
			 
			     A1=arco_S_P(xtrasl,ordinate_integrazione_S_2(1,ciclo))
           	     B1=arco_S_P(xtrasl,ordinate_integrazione_S_2(2,ciclo))
			 
	             alfa_j1=(B1-A1)/2.d0
	             beta_j1=(B1+A1)/2.d0

	             IF ((B1-A1).gt.10.d-15) THEN
	
	                 p2a = 0.d0

                     DO kj=1,Ngauss

	                     xinttrasl=(xint(kj)+1.d0)*0.5d0
	                     s=fi1(iplog,iqlog,xinttrasl)
	                     ds=dfi1(iplog,iqlog,xinttrasl)
	                     serv=alfa_j1*(2.d0*s-1.d0)+beta_j1
                         r2_1=(estremo_serv-serv)**2+xtrasl**2-2*(estremo_serv-serv)*xtrasl*COS(alpha)
						 r_Vucont(1)=(estremo_serv-serv)*COS(theta)+xtrasl*COS(theta2)
						 r_Vucont(2)=(estremo_serv-serv)*SIN(theta)+xtrasl*SIN(theta2)

	                     p2a = p2a + wint(kj)*ds*fiU(l_serv,serv,estremo_serv,grado_q)*((cp**2-cs**2)/(cp*cs))*delta_x*(cp*dsqrt((cs**2)*(delta_x**2)-r2_1)+cs*dsqrt((cp**2)*(delta_x**2)-r2_1))**(-1)*(r_Vucont(indice_i)*r_Vucont(indice_j)/r2_1-delta_kronecker(indice_i,indice_j)/2.d0)
	    
	                     IF (delta_kronecker(indice_i,indice_j).eq.1.d0) THEN
	    
	                         p2a= p2a+wint(kj)*ds*fiU(l_serv,serv,estremo_serv,grado_q)*coeff_3*((1.d0/(cp**2))*dlog(cp*delta_x+dsqrt(dabs(-r2_1+cp**2*delta_x**2)))+(1.d0/(cs**2))*dlog(cs*delta_x+dsqrt(dabs(-r2_1+cs**2*delta_x**2))))
	                     ENDIF
		  
                     END DO

        	     ENDIF

	             p2 = p2 +p2a*alfa_j1*w(ki)*fiU(l_var,xtrasl,estremo_var,grado_q)

             END DO
             int_2=int_2+p2*alfa
  
  
             !primo integrale     
             !-----------------

	         p3 = 0.d0
    
             IF (delta_kronecker(indice_i,indice_j).eq.1.d0) THEN
				 
                 DO ki=1,Ngauss
	     
	                 p3a = 0.d0

	                 xtrasl=alfa*x(ki)+beta	 
					 
				     A1=arco_S_P(xtrasl,ordinate_integrazione_S_2(1,ciclo))
           	         B1=arco_S_P(xtrasl,ordinate_integrazione_S_2(2,ciclo))

	                 alfa_j1=(B1-A1)/2.d0
	                 beta_j1=(B1+A1)/2.d0
	  
	                 if ((B1-A1).gt.10.d-15) then
					     
					     DO kj=1,Ngauss
						     iplog=1
							 iqlog=1
					         xinttrasl=(xint(kj)+1.d0)*0.5d0
	                         s=fi1(iplog,iqlog,xinttrasl)
	                         ds=dfi1(iplog,iqlog,xinttrasl)
	                         serv=alfa_j1*(2.d0*s-1.d0)+beta_j1
                             r2_1=(estremo_serv-serv)**2+xtrasl**2-2*(estremo_serv-serv)*xtrasl*COS(alpha)
						     
							 p3a=p3a+wint(kj)*ds*fiU(l_serv,serv,estremo_serv,grado_q)*dlog(dsqrt(r2_1))
							 
						 END DO
						 
	                 endif

                     p3=p3+p3a*alfa_j1*w(ki)*fiU(l_var,xtrasl,estremo_var,grado_q)
                 END DO
    
             ENDIF  
             int_2=int_2+coeff_2*p3*alfa
		 
	     ENDIF
     END DO	 
	 Vucont_sub_S=int_1-int_2
	 ! write(*,*) 'integrale in S', Vucont_sub_S 
	 
! !!!!!!!!!!!!!!!!!!!!!!ciclo integrazione su insieme ESP!!!!!!!!!!!!!!	 
     int_1=0.d0
	 int_2=0.d0
	 
     DO ciclo=1,10 !10 è la dimensione massima degli array ascisse.../ordinate...	
  	     if(ascisse_integrazione_SP_1(1,ciclo).ge.estremo_var)then
	         int_1=int_1+0.d0
         elseif(ascisse_integrazione_SP_1(1,ciclo) == ascisse_integrazione_SP_1(2,ciclo))then
	         int_1=int_1+0.d0
         else
	         if(ascisse_integrazione_SP_1(2,ciclo).ge.estremo_var)then
	             ascisse_integrazione_SP_1(2,ciclo)=estremo_var
	         endif
	
	         alfa=(ascisse_integrazione_SP_1(2,ciclo)-ascisse_integrazione_SP_1(1,ciclo))/2.d0
	         beta=(ascisse_integrazione_SP_1(2,ciclo)+ascisse_integrazione_SP_1(1,ciclo))/2.d0

             !secondo integrale   
             !------------------
	
             p2 = 0.d0
	         
			 if(ordinate_integrazione_SP_1(1,ciclo).eq.1)then
			     iplog=3
				 iqlog=1
			 endif
             if(ordinate_integrazione_SP_1(2,ciclo).eq.4)then
			     iplog=1
				 iqlog=3
			 endif
			 
	         DO ki=1,Ngauss 
	
                 xtrasl=alfa*x(ki)+beta
			 
			     A1=arco_S_P(xtrasl,ordinate_integrazione_SP_1(1,ciclo))
           	     B1=arco_S_P(xtrasl,ordinate_integrazione_SP_1(2,ciclo))
			 
	             alfa_j1=(B1-A1)/2.d0
	             beta_j1=(B1+A1)/2.d0
                 
				 p2a = 0.d0
				 
                 IF ((B1-A1).gt.10.d-15) THEN
	
                     DO kj=1,Ngauss

	                     xinttrasl=(xint(kj)+1.d0)*0.5d0
	                     s=fi1(iplog,iqlog,xinttrasl)
	                     ds=dfi1(iplog,iqlog,xinttrasl)
	                     serv=alfa_j1*(2.d0*s-1.d0)+beta_j1
                         r2_1=(estremo_serv-serv)**2+xtrasl**2-2*(estremo_serv-serv)*xtrasl*COS(alpha)
						 r_Vucont(1)=(estremo_serv-serv)*COS(theta)+xtrasl*COS(theta2)
						 r_Vucont(2)=(estremo_serv-serv)*SIN(theta)+xtrasl*SIN(theta2)

	                     p2a = p2a + wint(kj)*ds*fiU(l_serv,serv,estremo_serv,grado_q)*((cp**2*r2_1)**(-1))*cp*delta_x*dsqrt(dabs(-r2_1+cp**2*delta_x**2))*(r_Vucont(indice_i)*r_Vucont(indice_j)/r2_1-delta_kronecker(indice_i,indice_j)/2.d0)
	    
	                     IF (delta_kronecker(indice_i,indice_j).eq.1.d0) THEN
	    
	                         p2a= p2a+wint(kj)*ds*fiU(l_serv,serv,estremo_serv,grado_q)*coeff_5*(1.d0/(cp**2))*dlog(cp*delta_x+dsqrt(dabs(-r2_1+cp**2*delta_x**2)))
	                     ENDIF						 

                     END DO

        	     ENDIF

	             p2 = p2 +p2a*alfa_j1*w(ki)*fiU(l_var,xtrasl,estremo_var,grado_q)

             END DO
             int_1=int_1+p2*alfa
  
  
             !primo integrale     
             !-----------------

	         p3 = 0.d0
    
             IF (delta_kronecker(indice_i,indice_j).eq.1.d0) THEN
               
                 DO ki=1,Ngauss
	     
	                 p3a = 0.d0
	 
	                 xtrasl=alfa*x(ki)+beta
	             
				     A1=arco_S_P(xtrasl,ordinate_integrazione_SP_1(1,ciclo))
           	         B1=arco_S_P(xtrasl,ordinate_integrazione_SP_1(2,ciclo))

	                 alfa_j1=(B1-A1)/2.d0
	                 beta_j1=(B1+A1)/2.d0
	  
	                 if ((B1-A1).gt.10.d-15) then
					     
					     DO kj=1,Ngauss
						     iplog=1
							 iqlog=1
					         xinttrasl=(xint(kj)+1.d0)*0.5d0
	                         s=fi1(iplog,iqlog,xinttrasl)
	                         ds=dfi1(iplog,iqlog,xinttrasl)
	                         serv=alfa_j1*(2.d0*s-1.d0)+beta_j1
                             r2_1=(estremo_serv-serv)**2+xtrasl**2-2*(estremo_serv-serv)*xtrasl*COS(alpha)
						     
							 p3a=p3a+wint(kj)*ds*fiU(l_serv,serv,estremo_serv,grado_q)*dlog(dsqrt(r2_1))
							 
						 END DO
						 
	                 endif

                     p3=p3+p3a*alfa_j1*w(ki)*fiU(l_var,xtrasl,estremo_var,grado_q)
                 END DO
    
             ENDIF  
			 
			 int_1=int_1-(1.d0/(cp**2))*coeff_5*p3*alfa
		 
	     ENDIF
     END DO  

     DO ciclo=1,10 !10 è la dimensione massima degli array ascisse.../ordinate...	
  	     if(ascisse_integrazione_SP_2(1,ciclo).ge.estremo_var)then
	         int_2=int_2+0.d0
         elseif(ascisse_integrazione_SP_2(1,ciclo) == ascisse_integrazione_SP_2(2,ciclo))then
	         int_2=int_2+0.d0
         else
	         if(ascisse_integrazione_SP_2(2,ciclo).ge.estremo_var)then
	             ascisse_integrazione_SP_2(2,ciclo)=estremo_var
	         endif
	
	         alfa=(ascisse_integrazione_SP_2(2,ciclo)-ascisse_integrazione_SP_2(1,ciclo))/2.d0
	         beta=(ascisse_integrazione_SP_2(2,ciclo)+ascisse_integrazione_SP_2(1,ciclo))/2.d0

             !secondo integrale   
             !------------------
	
             p2 = 0.d0
	         
			 if(ordinate_integrazione_SP_2(1,ciclo).eq.1)then
			     iplog=3
				 iqlog=1
			 endif
             if(ordinate_integrazione_SP_2(2,ciclo).eq.4)then
			     iplog=1
				 iqlog=3
			 endif
			 
	         DO ki=1,Ngauss 
	
                 xtrasl=alfa*x(ki)+beta
			 
			     A1=arco_S_P(xtrasl,ordinate_integrazione_SP_2(1,ciclo))
           	     B1=arco_S_P(xtrasl,ordinate_integrazione_SP_2(2,ciclo))
			 
	             alfa_j1=(B1-A1)/2.d0
	             beta_j1=(B1+A1)/2.d0

	             IF ((B1-A1).gt.10.d-15) THEN
	
	                 p2a = 0.d0

                     DO kj=1,Ngauss

	                     xinttrasl=(xint(kj)+1.d0)*0.5d0
	                     s=fi1(iplog,iqlog,xinttrasl)
	                     ds=dfi1(iplog,iqlog,xinttrasl)
	                     serv=alfa_j1*(2.d0*s-1.d0)+beta_j1
                         r2_1=(estremo_serv-serv)**2+xtrasl**2-2*(estremo_serv-serv)*xtrasl*COS(alpha)
						 r_Vucont(1)=(estremo_serv-serv)*COS(theta)+xtrasl*COS(theta2)
						 r_Vucont(2)=(estremo_serv-serv)*SIN(theta)+xtrasl*SIN(theta2)
						 
	                     p2a = p2a + wint(kj)*ds*fiU(l_serv,serv,estremo_serv,grado_q)*((cp**2*r2_1)**(-1))*cp*delta_x*dsqrt(dabs(-r2_1+cp**2*delta_x**2))*(r_Vucont(indice_i)*r_Vucont(indice_j)/r2_1-delta_kronecker(indice_i,indice_j)/2.d0)
	    
	                     IF (delta_kronecker(indice_i,indice_j).eq.1.d0) THEN
	    
	                         p2a= p2a+wint(kj)*ds*fiU(l_serv,serv,estremo_serv,grado_q)*coeff_5*(1.d0/(cp**2))*dlog(cp*delta_x+dsqrt(dabs(-r2_1+cp**2*delta_x**2)))
	                     ENDIF							 
		  
                     END DO

        	     ENDIF

	             p2 = p2 +p2a*alfa_j1*w(ki)*fiU(l_var,xtrasl,estremo_var,grado_q)

             END DO
             int_2=int_2+p2*alfa
  
  
             !primo integrale     
             !-----------------

	         p3 = 0.d0
    
             IF (delta_kronecker(indice_i,indice_j).eq.1.d0) THEN
    
                 DO ki=1,Ngauss
	     
	                 p3a = 0.d0
	 
	                 xtrasl=alfa*x(ki)+beta
	             
				     A1=arco_S_P(xtrasl,ordinate_integrazione_SP_2(1,ciclo))
           	         B1=arco_S_P(xtrasl,ordinate_integrazione_SP_2(2,ciclo))

	                 alfa_j1=(B1-A1)/2.d0
	                 beta_j1=(B1+A1)/2.d0
	  
	                 if ((B1-A1).gt.10.d-15) then
					     
					     DO kj=1,Ngauss
						     iplog=1
							 iqlog=1
					         xinttrasl=(xint(kj)+1.d0)*0.5d0
	                         s=fi1(iplog,iqlog,xinttrasl)
	                         ds=dfi1(iplog,iqlog,xinttrasl)
	                         serv=alfa_j1*(2.d0*s-1.d0)+beta_j1
                             r2_1=(estremo_serv-serv)**2+xtrasl**2-2*(estremo_serv-serv)*xtrasl*COS(alpha)
						     
							 p3a=p3a+wint(kj)*ds*fiU(l_serv,serv,estremo_serv,grado_q)*dlog(dsqrt(r2_1))
							 
						 END DO
						 
	                 endif

                     p3=p3+p3a*alfa_j1*w(ki)*fiU(l_var,xtrasl,estremo_var,grado_q)
                 END DO
    
             ENDIF  
	         int_2=int_2-(1.d0/(cp**2))*coeff_5*p3*alfa	 
		 
	     ENDIF
     END DO	 
	 Vucont_sub_P=int_1-int_2
     ! write(*,*) 'integrale in S', Vucont_sub_P
	 ! pause
 endif
     Vucont_sub=Vucont_sub_S+Vucont_sub_P
	  
  RETURN

   CONTAINS 
   
   FUNCTION retta_S_P(x_linea,tipo_segmento) !!! per definire le rette dei domini in caso di elementi allineati (o eventualmete limiti delle ordinate)
      REAL(kind=8) :: retta_S_P, x_linea
      INTEGER(kind=4) :: tipo_segmento
	  select case (tipo_segmento)
         case (0) 
             retta_S_P=0.d0 
         case (1)
             retta_S_P = x_linea+estremo_serv-cp*delta_x
         case (2) 
             retta_S_P = x_linea+estremo_serv-cs*delta_x 
         case (3)
             retta_S_P=estremo_serv
         case default
         print*, "entrata non valida" 
      end select
       
   END FUNCTION retta_S_P
   
   FUNCTION arco_S_P(x_linea,tipo_arco)
     REAL(kind=8) :: arco_S_P, x_linea
	 INTEGER(kind=4) :: tipo_arco
	 select case (tipo_arco)
	     case(0)
		     arco_S_P=0.d0
		 case(1)
		     arco_S_P=estremo_serv-x_linea*COS(alpha)-dsqrt(-(x_linea*SIN(alpha))**2+cp**2*delta_x**2)
		 case(2)
		     arco_S_P=estremo_serv-x_linea*COS(alpha)-dsqrt(-(x_linea*SIN(alpha))**2+cs**2*delta_x**2)
		 case(3)
		     arco_S_P=estremo_serv-x_linea*COS(alpha)+dsqrt(-(x_linea*SIN(alpha))**2+cs**2*delta_x**2)
		 case(4)
             arco_S_P=estremo_serv-x_linea*COS(alpha)+dsqrt(-(x_linea*SIN(alpha))**2+cp**2*delta_x**2)		 
		 case(5)
		     arco_S_P=estremo_serv
		 case default
		     print*, "entrata non valida" 
	 end select
   END FUNCTION arco_S_P
  END
    