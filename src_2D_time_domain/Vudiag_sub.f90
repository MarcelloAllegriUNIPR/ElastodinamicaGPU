      Double precision function Vudiag_sub(li,lj,i,tempo1,tempo2,indice_i,indice_j)
!
!     INTEGRALE DELLE FUNZIONI DI FORMA li,lj SULL'ELEMENTO i
!     ELEMENTI SOVRAPPOSTI
!
  USE variable_2Dgeneral
  
  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABILI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  INTEGER(kind=4),INTENT(IN):: i, li, lj, indice_i, indice_j
  REAL(kind=8),INTENT(IN):: tempo1, tempo2
  
  !Variabili locali
  INTEGER(kind=4):: ind_gauss, Ngauss, AllocateStatus, iplog, iqlog, ki, kj, ip, iq, Ngaussi
  REAL(kind=8):: Vudiag_sub11, Vudiag_sub12, Vudiag_sub21, Vudiag_sub22, Vudiag_sub31, Vudiag_sub32, delta_x, xi1, xi2, xm1, xm2, xm, xn, alfa, alfa2, beta, & 
    p2, A1, B1, alfa_j1, beta_j1, p2a, s, ds, xtrasl, xinttrasl, r2_1, p3, p3a, cs1, serv, sm1, sm2, sn, sm, Vudiag2_sub11, Vudiag2_sub12, Vudiag2_sub21, Vudiag2_sub22, &
    Vudiag2_sub31, Vudiag2_sub32, Vudiag2_sub41, Vudiag2_sub42
  REAL(kind=8):: estremo, th, tk, theta, coeff_1, coeff_2, coeff_3, coeff_4, coeff_5, cp, cs
  INTEGER(kind=4), DIMENSION(2,2) :: delta_kronecker
  REAL(kind=8), DIMENSION(2) :: sigma
  REAL(kind=8),DIMENSION(:),ALLOCATABLE:: x, w, xint, wint, ti, vi, ww1
  REAL(kind=8),EXTERNAL:: fi1, dfi1, fiU
    
  !!!!!!!!!!!!!!!!!!!!!!!!!! CORPO della SUBROUTINE !!!!!!!!!!!!!!!!!!
  cs=velC_S
  cp=velC_P
  
  estremo=list_elements(i)%length
  theta=ACOS((list_nodes(i+1)%coordinates(1)-list_nodes(i)%coordinates(1))/estremo)
  if(list_nodes(i+1)%coordinates(2).lt.list_nodes(i)%coordinates(2))then
     theta=-theta
  endif
  
  
  delta_kronecker(1,1)=1
  delta_kronecker(1,2)=0
  delta_kronecker(2,1)=0
  delta_kronecker(2,2)=1
  sigma(1)=COS(theta)
  sigma(2)=SIN(theta)
  
  coeff_1=((cp**2-cs**2)/(cp*cs))*(sigma(indice_i)*sigma(indice_j)-delta_kronecker(indice_i,indice_j)/2.d0)
  coeff_2=-((cp**2+cs**2)/(cp**2*cs**2))*delta_kronecker(indice_i,indice_j)/2.d0
  coeff_3=delta_kronecker(indice_i,indice_j)/2.d0
  coeff_4=sigma(indice_i)*sigma(indice_j)-delta_kronecker(indice_i,indice_j)/2.d0
  coeff_5=delta_kronecker(indice_i,indice_j)/2.d0
	
  Vudiag_sub11=0.d0
  Vudiag_sub12=0.d0
  Vudiag_sub21=0.d0
  Vudiag_sub22=0.d0
  Vudiag_sub31=0.d0
  Vudiag_sub32=0.d0
  
  Vudiag2_sub11=0.d0
  Vudiag2_sub12=0.d0
  Vudiag2_sub21=0.d0
  Vudiag2_sub22=0.d0
  Vudiag2_sub31=0.d0
  Vudiag2_sub32=0.d0
  Vudiag2_sub41=0.d0
  Vudiag2_sub42=0.d0
  
  Vudiag_sub=0.d0
	
  delta_x=tempo1-tempo2
  IF (delta_x.le.0.d0) RETURN
  
  if(cs*delta_x.le.estremo)then
  ind_gauss=6
  else
  ind_gauss=6
  endif
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

  iplog=2
  iqlog=2
	
  xi1=0.d0
  xi2=estremo
												
  xm1=estremo-cs*delta_x
  xm2=cs*delta_x

  if(cs*delta_x.ge.xi2)then
  xm=xi1
  xn=xi1
  elseif(cs*delta_x.lt.xi2/2.d0)then
  xm=xm2
  xn=xm1
  elseif(cs*delta_x.eq.xi2/2.d0)then
  xm=xm2
  xn=xm
  elseif(cs*delta_x.gt.xi2/2.d0)then
  xm=xm1
  xn=xm2
  endif
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
    
	! p2=0.d0
  	! DO ki=1,Ngauss 
	
      ! xtrasl=estremo*(x(ki)+1.d0)/2.d0

	  ! p2a = 0.d0

        ! DO kj=1,Ngauss

	        ! xinttrasl=estremo*(xint(kj)+1.d0)/2.d0

	        ! p2a = p2a + wint(kj)*fiU(lj,xinttrasl,estremo,grado_q)
		  
        ! END DO

	  ! p2 = p2 +p2a*w(ki)*fiU(li,xtrasl,estremo,grado_q)

    ! END DO
	! Vudiag_sub=p2*estremo/4.d0
  
  
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! write(*,*) 'angolo theta', theta
  ! write(*,*) 'indici i e j', indice_i, indice_j
  ! write(*,*) 'elemento', i

!______________________________________
!	Sul primo sottointervallo
	
	if(xm-xi1.gt.0)then
	alfa2=xm-xi1
	alfa=alfa2/2.d0
	beta=(xm+xi1)/2.d0

!  secondo integrale   
!------------------
	iplog=1
    iqlog=2
	
	p2 = 0.d0
	 
	DO ki=1,Ngauss 
	
      xtrasl=alfa*(x(ki)+1.d0)
	  A1=dmax1(0.d0,-cs*delta_x+xtrasl+xi1)
	  B1=dmin1(xi2,cs*delta_x+xtrasl+xi1)

	  alfa_j1=(B1-A1)/2.d0
	  beta_j1=(B1+A1)/2.d0

	  p2a = 0.d0
	  
	  IF ((B1-A1).gt.10.d-15) THEN

        DO kj=1,Ngauss

	        xinttrasl=(xint(kj)+1.d0)*0.5d0
	        s=fi1(iplog,iqlog,xinttrasl)
	        ds=dfi1(iplog,iqlog,xinttrasl)
	        serv=alfa_j1*(2.d0*s-1.d0)+beta_j1

		    r2_1=(xtrasl-serv+xi1)**2

	        p2a = p2a + wint(kj)*ds*coeff_1*delta_x*(cp*dsqrt((cs**2)*(delta_x**2)-r2_1)+cs*dsqrt((cp**2)*(delta_x**2)-r2_1))**(-1)*&
			                                 fiU(lj,serv,estremo,grado_q)*fiU(li,alfa*x(ki)+beta,estremo,grado_q)
	        
	         IF (delta_kronecker(indice_i,indice_j).eq.1.d0) THEN
	    
	         p2a= p2a+wint(kj)*ds*coeff_3*((1.d0/(cp**2))*dlog(cp*delta_x+dsqrt(dabs(-r2_1+cp**2*delta_x**2)))+&
			             (1.d0/(cs**2))*dlog(cs*delta_x+dsqrt(dabs(-r2_1+cs**2*delta_x**2))))*fiU(lj,serv,estremo,grado_q)*fiU(li,alfa*x(ki)+beta,estremo,grado_q)
	         ENDIF
		  
        END DO

	  ENDIF

	  p2 = p2 +p2a*alfa_j1*w(ki)

    END DO
   Vudiag_sub11= p2*alfa
  
!   primo integrale     
!   -----------------
	 
	ip=2
	iq=1

	  p3 = 0.d0
    
    IF (delta_kronecker(indice_i,indice_j).eq.1.d0) THEN
    
    DO ki=1,Ngauss
	    
	    p3a = 0.d0
	 
	    xtrasl=0.5d0*x(ki)+0.5d0
	    serv=alfa2*fi1(ip,iq,xtrasl)+xi1
	    A1=dmax1(0.d0,-cs*delta_x+serv)
	    B1=dmin1(xi2,cs*delta_x+serv)

	    alfa_j1=(B1-A1)/2.d0
	    beta_j1=(B1+A1)/2.d0
	  
	    if ((B1-A1).gt.10.d-15) then

		    cs1=(serv-beta_j1)/alfa_j1

	        if (dabs(cs1).lt.1.25d0) then

		        call mopelog(Ngaussi,ti,cs1,0.d0,ww1,1)
    	 
                DO kj=1,Ngaussi
			        p3a=p3a+vi(kj)*ww1(kj)*fiU(lj,alfa_j1*ti(kj)+beta_j1,estremo,grado_q)*fiU(li,serv,estremo,grado_q)
                END DO

	            DO kj=1,Ngauss
	                p3a=p3a+w(kj)*dlog(alfa_j1)*fiU(lj,alfa_j1*x(kj)+beta_j1,estremo,grado_q)*fiU(li,serv,estremo,grado_q)
                END DO

	        else
    	        Write(*,*) 'non ho usato la mopelog 1'
				pause
	            DO kj=1,Ngauss
	                p3a=p3a+w(kj)*(dlog(dabs(cs1-x(kj)))+dlog(alfa_j1))*fiU(lj,alfa_j1*x(kj)+beta_j1,estremo,grado_q)*fiU(li,serv,estremo,grado_q)
                END DO
                
	        endif

	    endif

        p3=p3+p3a*alfa_j1*w(ki)*dfi1(ip,iq,xtrasl)!

    END DO
    
     ENDIF  
  
	  Vudiag_sub=Vudiag_sub+(p2+coeff_2*p3)*alfa
     Vudiag_sub12=coeff_2*p3*alfa
  	  endif
! !________________________________________________
!	Sul secondo sottointervallo

	if(xn-xm.gt.0)then
	alfa2=xn-xm
	alfa=alfa2/2.d0
	beta=(xm+xn)/2.d0

!   secondo integrale   
! ------------------
	iplog=2
    iqlog=2
	
	p2 = 0.d0
	 
	DO ki=1,Ngauss 
	
	  xtrasl=alfa*(x(ki)+1.d0)
	  A1=dmax1(0.d0,-cs*delta_x+xtrasl+xm)
	  B1=dmin1(xi2,cs*delta_x+xtrasl+xm)

	  alfa_j1=(B1-A1)/2.d0
	  beta_j1=(B1+A1)/2.d0

	  p2a = 0.d0
	  
	    if ((B1-A1).gt.10.d-15) then

            DO kj=1,Ngauss

	            xinttrasl=(xint(kj)+1.d0)*0.5d0
	            s=fi1(iplog,iqlog,xinttrasl)
	            ds=dfi1(iplog,iqlog,xinttrasl)
	            serv=alfa_j1*(2.d0*s-1.d0)+beta_j1

		        r2_1=(xtrasl-serv+xm)**2

	            p2a = p2a + wint(kj)*ds*fiU(lj,serv,estremo,grado_q)*fiU(li,alfa*x(ki)+beta,estremo,grado_q)*coeff_1*delta_x&
				                       *(cp*dsqrt((cs**2)*(delta_x**2)-r2_1)+cs*dsqrt((cp**2)*(delta_x**2)-r2_1))**(-1)
	    
	            IF (delta_kronecker(indice_i,indice_j).eq.1.d0) THEN
	    
	            p2a= p2a+wint(kj)*ds*fiU(lj,serv,estremo,grado_q)*fiU(li,alfa*x(ki)+beta,estremo,grado_q)*coeff_3*((1.d0/(cp**2))&
				           *dlog(cp*delta_x+dsqrt(dabs(-r2_1+cp**2*delta_x**2)))+(1.d0/(cs**2))*dlog(cs*delta_x+dsqrt(dabs(-r2_1+cs**2*delta_x**2))))
	    
    	        ENDIF
		  
            END DO

	    endif

	    p2 = p2 +p2a*alfa_j1*w(ki)

    END DO
	Vudiag_sub21=p2*alfa
  
!   primo integrale     
!   -----------------
	 
	 ip=2
	 iq=2

	  p3 = 0.d0
	  
	 IF (delta_kronecker(indice_i,indice_j).eq.1.d0) THEN
	     
	 DO ki=1,Ngauss
	    
	    p3a = 0.d0
	 
	    xtrasl=0.5d0*x(ki)+0.5d0
	    serv=alfa2*fi1(ip,iq,xtrasl)+xm
        A1=dmax1(0.d0,-cs*delta_x+serv)
	    B1=dmin1(xi2,cs*delta_x+serv)

	    alfa_j1=(B1-A1)/2.d0
	    beta_j1=(B1+A1)/2.d0
	  
	    if ((B1-A1).gt.10.d-15) then

		    cs1=(serv-beta_j1)/alfa_j1

	        if (dabs(cs1).lt.1.25d0) then

		        call mopelog(Ngaussi,ti,cs1,0.d0,ww1,1)
	 
                DO kj=1,Ngaussi
		            p3a=p3a+vi(kj)*ww1(kj)*fiU(lj,alfa_j1*ti(kj)+beta_j1,estremo,grado_q)*fiU(li,serv,estremo,grado_q)
                END DO

	            DO kj=1,Ngauss
	                p3a=p3a+w(kj)*dlog(alfa_j1)*fiU(lj,alfa_j1*x(kj)+beta_j1,estremo,grado_q)*fiU(li,serv,estremo,grado_q)
                END DO

	        else
	            Write(*,*) 'non ho usato la mopelog 2'
				pause
	            DO kj=1,Ngauss
	                p3a=p3a+w(kj)*(dlog(dabs(cs1-x(kj)))+dlog(alfa_j1))*fiU(lj,alfa_j1*x(kj)+beta_j1,estremo,grado_q)*fiU(li,serv,estremo,grado_q)
                END	DO
                
            endif

	    endif

        p3=p3+p3a*alfa_j1*w(ki)*dfi1(ip,iq,xtrasl)!

    END DO

    ENDIF 

     Vudiag_sub=Vudiag_sub+(p2+coeff_2*p3)*alfa

    Vudiag_sub22=coeff_2*p3*alfa

    
	 endif
	
!________________________________________________
!	Sul terzo sottointervallo

	if(xi2-xn.gt.0)then
	alfa2=xi2-xn
	alfa=alfa2/2.d0
	beta=(xi2+xn)/2.d0

!   secondo integrale   
! ------------------
	iplog=2
    iqlog=1
	
	if(xn.eq.xi1)then
	 iqlog=2
	endif
	
	p2 = 0.d0
	 
	DO ki=1,Ngauss 

	    xtrasl=alfa*(x(ki)+1.d0)
        A1=dmax1(0.d0,-cs*delta_x+xtrasl+xn)
	    B1=dmin1(xi2,cs*delta_x+xtrasl+xn)

	    alfa_j1=(B1-A1)/2.d0
	    beta_j1=(B1+A1)/2.d0

		p2a = 0.d0
		
	    if ((B1-A1).gt.10.d-15) then

            DO kj=1,Ngauss

	            xinttrasl=(xint(kj)+1.d0)*0.5d0
	            s=fi1(iplog,iqlog,xinttrasl)
	            ds=dfi1(iplog,iqlog,xinttrasl)
	            serv=alfa_j1*(2.d0*s-1.d0)+beta_j1
	            
                r2_1=(xtrasl-alfa_j1*(2.d0*s-1.d0)-beta_j1+xn)**2

	            p2a = p2a + wint(kj)*ds*fiU(lj,serv,estremo,grado_q)*fiU(li,alfa*x(ki)+beta,estremo,grado_q)*coeff_1*delta_x&
				                      *(cp*dsqrt((cs**2)*(delta_x**2)-r2_1)+cs*dsqrt((cp**2)*(delta_x**2)-r2_1))**(-1)
	    
	            IF (delta_kronecker(indice_i,indice_j).eq.1.d0) THEN
	    
	            p2a= p2a + wint(kj)*ds*fiU(lj,serv,estremo,grado_q)*fiU(li,alfa*x(ki)+beta,estremo,grado_q)*coeff_3*((1.d0/(cp**2))&
				       *dlog(cp*delta_x+dsqrt(dabs(-r2_1+cp**2*delta_x**2)))+(1.d0/(cs**2))*dlog(cs*delta_x+dsqrt(dabs(-r2_1+cs**2*delta_x**2))))
	    
	            ENDIF
		  
            END DO

	    endif

	    p2 = p2 +p2a*alfa_j1*w(ki)

    END DO
    
    Vudiag_sub31=p2*alfa  		
 
!   primo integrale     
!   -----------------
	 
	 ip=1
	 iq=2
	 
	 if(xn.eq.xi1)then
	 ip=2
	 endif

	  p3 = 0.d0
	 
	 IF (delta_kronecker(indice_i,indice_j).eq.1.d0) THEN
	 
	 DO ki=1,Ngauss
	 
	    p3a = 0.d0
	 
	    xtrasl=0.5d0*x(ki)+0.5d0
	    serv=alfa2*fi1(ip,iq,xtrasl)+xn
        A1=dmax1(0.d0,-cs*delta_x+serv)
	    B1=dmin1(xi2,cs*delta_x+serv)

	    alfa_j1=(B1-A1)/2.d0
	    beta_j1=(B1+A1)/2.d0
	  
	    if ((B1-A1).gt.10.d-15) then
	        
	        cs1=(serv-beta_j1)/alfa_j1

	        if (dabs(cs1).lt.1.25d0) then

		        call mopelog(Ngaussi,ti,cs1,0.d0,ww1,1)
	 
                DO kj=1,Ngaussi
                    p3a=p3a+vi(kj)*ww1(kj)*fiU(lj,alfa_j1*ti(kj)+beta_j1,estremo,grado_q)*fiU(li,serv,estremo,grado_q)
                END DO

                DO kj=1,Ngauss
	                p3a=p3a+w(kj)*dlog(alfa_j1)*fiU(lj,alfa_j1*x(kj)+beta_j1,estremo,grado_q)*fiU(li,serv,estremo,grado_q)
                END DO

	        else
	            Write(*,*) 'non ho usato la mopelog 3'
				pause
	            DO kj=1,Ngauss
	                p3a=p3a+w(kj)*(dlog(dabs(cs1-x(kj)))+dlog(alfa_j1))*fiU(lj,alfa_j1*x(kj)+beta_j1,estremo,grado_q)*fiU(li,serv,estremo,grado_q)
	            END DO
	            
	        endif

	    endif

        p3=p3+p3a*alfa_j1*w(ki)*dfi1(ip,iq,xtrasl)!

    END DO
    
    ENDIF 

	 Vudiag_sub=Vudiag_sub+(p2+coeff_2*p3)*alfa
    Vudiag_sub32=coeff_2*p3*alfa
	 endif  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!DOMINI AGGIUNTIVI!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  												
  sm1=estremo-cp*delta_x
  sm2=cp*delta_x
  
  if(sm1.gt.0.d0)then
  sm=sm1
  sn=sm2
  elseif(sm1.eq.0.d0)then
  sm=sm1
  sn=sm2
  elseif(sm1.lt.0.d0)then
  sm=0
  sn=estremo
  endif
  
  !xm1=estremo-cs*delta_x
  !xm2=cs*delta_x
  
  if(xm2.gt.estremo)then
  xm=0
  xn=estremo
  elseif(xm2.eq.estremo)then
  xm=0
  xn=estremo
  elseif(xm2.lt.estremo)then
  xm=xm1
  xn=xm2
  endif
  
  
  !______________________________________
!	Su intervalli 1 e 4
	
	if(sm.gt.0)then
!____________________________________intervallo 1__________________________

	alfa=sm/2.d0
	beta=sm/2.d0

!   secondo integrale   
! ------------------
	
	iplog=1
    iqlog=2
	
	p2 = 0.d0
	 
	DO ki=1,Ngauss 
	
      xtrasl=alfa*x(ki)+beta
	  A1=cs*delta_x+xtrasl
	  B1=cp*delta_x+xtrasl

	  alfa_j1=(B1-A1)/2.d0
	  beta_j1=(B1+A1)/2.d0

	  p2a = 0.d0
	  
	  IF ((B1-A1).gt.10.d-15) THEN

        DO kj=1,Ngauss
            
            
            xinttrasl=(xint(kj)+1.d0)*0.5d0
	        s=fi1(iplog,iqlog,xinttrasl)
	        ds=dfi1(iplog,iqlog,xinttrasl)
	        serv=alfa_j1*(2.d0*s-1.d0)+beta_j1
            r2_1=(xtrasl-serv)**2

	        p2a = p2a + wint(kj)*ds*fiU(lj,serv,estremo,grado_q)*coeff_4*((cp**2*r2_1)**(-1))*cp*delta_x*dsqrt(dabs(-r2_1+cp**2*delta_x**2))
	    
	        IF (delta_kronecker(indice_i,indice_j).eq.1.d0) THEN
	    
	        p2a= p2a+wint(kj)*ds*fiU(lj,serv,estremo,grado_q)*coeff_5*(1.d0/(cp**2))*dlog(cp*delta_x+dsqrt(dabs(-r2_1+cp**2*delta_x**2)))
	        ENDIF
		  
        END DO

	  ENDIF

	  p2 = p2 +p2a*alfa_j1*w(ki)*fiU(li,xtrasl,estremo,grado_q)

    END DO
    Vudiag2_sub11= p2*alfa
    Vudiag_sub=Vudiag_sub+Vudiag2_sub11
    
  
!   primo integrale     
!   -----------------

	p3 = 0.d0
	
    IF (delta_kronecker(indice_i,indice_j).eq.1.d0) THEN
    
    DO ki=1,Ngauss
	    
	    p3a = 0.d0
	 
	    xtrasl=alfa*x(ki)+beta
	    A1=cs*delta_x+xtrasl
	    B1=cp*delta_x+xtrasl

	    alfa_j1=(B1-A1)/2.d0
	    beta_j1=(B1+A1)/2.d0
	  
	    if ((B1-A1).gt.10.d-15) then

		    cs1=(xtrasl-beta_j1)/alfa_j1

	        if (dabs(cs1).lt.1.25d0) then
                Write(*,*) 'ho usato la mopelog 1 cp'
				pause
				
		        call mopelog(Ngaussi,ti,cs1,0.d0,ww1,1)
    	 
                DO kj=1,Ngaussi
			        p3a=p3a+vi(kj)*ww1(kj)*fiU(lj,alfa_j1*ti(kj)+beta_j1,estremo,grado_q)
                END DO

	            DO kj=1,Ngauss
	                p3a=p3a+w(kj)*fiU(lj,alfa_j1*x(kj)+beta_j1,estremo,grado_q)*dlog(alfa_j1)
                END DO

	        else
    	        
	            DO kj=1,Ngauss
	                p3a=p3a+w(kj)*fiU(lj,alfa_j1*x(kj)+beta_j1,estremo,grado_q)*(dlog(dabs(cs1-x(kj)))+dlog(alfa_j1))
                END DO
                
	        endif

	    endif

        p3=p3+p3a*alfa_j1*w(ki)*fiU(li,xtrasl,estremo,grado_q)
    END DO
    
    ENDIF  
  
    Vudiag2_sub12=-(1.d0/(cp**2))*coeff_5*p3*alfa
    Vudiag_sub=Vudiag_sub+Vudiag2_sub12
	
  	
  	!____________________________________intervallo 4__________________________
  	
  	
  	alfa=(estremo-sn)/2.d0
	beta=(estremo+sn)/2.d0

!   secondo integrale   
! ------------------
	
	p2 = 0.d0
	 
	iplog=2
    iqlog=1
	 
	DO ki=1,Ngauss 
	
      xtrasl=alfa*x(ki)+beta
	  A1=-cp*delta_x+xtrasl
	  B1=-cs*delta_x+xtrasl

	  alfa_j1=(B1-A1)/2.d0
	  beta_j1=(B1+A1)/2.d0

	  IF ((B1-A1).gt.10.d-15) THEN
	
	    p2a = 0.d0

        DO kj=1,Ngauss

	        xinttrasl=(xint(kj)+1.d0)*0.5d0
	        s=fi1(iplog,iqlog,xinttrasl)
	        ds=dfi1(iplog,iqlog,xinttrasl)
	        serv=alfa_j1*(2.d0*s-1.d0)+beta_j1
            r2_1=(xtrasl-serv)**2

	        p2a = p2a + wint(kj)*ds*fiU(lj,serv,estremo,grado_q)*coeff_4*((cp**2*r2_1)**(-1))*cp*delta_x*dsqrt(dabs(-r2_1+cp**2*delta_x**2))
	    
	        IF (delta_kronecker(indice_i,indice_j).eq.1.d0) THEN
	    
	        p2a= p2a+wint(kj)*ds*fiU(lj,serv,estremo,grado_q)*coeff_5*(1.d0/(cp**2))*dlog(cp*delta_x+dsqrt(dabs(-r2_1+cp**2*delta_x**2)))
	        ENDIF
		  
        END DO

	  ENDIF

	  p2 = p2 +p2a*alfa_j1*w(ki)*fiU(li,xtrasl,estremo,grado_q)

    END DO
    Vudiag2_sub41= p2*alfa
    Vudiag_sub=Vudiag_sub+Vudiag2_sub41
    
  
!   primo integrale     
!   -----------------

	p3 = 0.d0
    
    IF (delta_kronecker(indice_i,indice_j).eq.1.d0) THEN
    
    DO ki=1,Ngauss
	    
	    p3a = 0.d0
	 
	    xtrasl=alfa*x(ki)+beta
	    A1=-cp*delta_x+xtrasl
	    B1=-cs*delta_x+xtrasl

	    alfa_j1=(B1-A1)/2.d0
	    beta_j1=(B1+A1)/2.d0
	  
	    if ((B1-A1).gt.10.d-15) then

		    cs1=(xtrasl-beta_j1)/alfa_j1

	        if (dabs(cs1).lt.1.25d0) then
                Write(*,*) 'ho usato la mopelog 4 cp'
				pause
				
		        call mopelog(Ngaussi,ti,cs1,0.d0,ww1,1)
    	 
                DO kj=1,Ngaussi
			        p3a=p3a+vi(kj)*ww1(kj)*fiU(lj,alfa_j1*ti(kj)+beta_j1,estremo,grado_q)
                END DO

	            DO kj=1,Ngauss
	                p3a=p3a+w(kj)*fiU(lj,alfa_j1*x(kj)+beta_j1,estremo,grado_q)*dlog(alfa_j1)
                END DO

	        else
    	        
	            DO kj=1,Ngauss
	                p3a=p3a+w(kj)*fiU(lj,alfa_j1*x(kj)+beta_j1,estremo,grado_q)*(dlog(dabs(cs1-x(kj)))+dlog(alfa_j1))
                END DO
                
	        endif

	    endif

        p3=p3+p3a*alfa_j1*w(ki)*fiU(li,xtrasl,estremo,grado_q)
    END DO
    
    ENDIF  
  
    Vudiag2_sub42=-(1.d0/(cp**2))*coeff_5*p3*alfa
    Vudiag_sub=Vudiag_sub+Vudiag2_sub42

  endif
  
  !______________________________________
!	Su intervalli 2 e 3
	
	if(sm.lt.xm)then
!____________________________________intervallo 2__________________________
	alfa=(xm-sm)/2.d0
	beta=(xm+sm)/2.d0
    
	iplog=1
    iqlog=2
	
!   secondo integrale   
! ------------------
	
	p2 = 0.d0
	 
	DO ki=1,Ngauss 
	
      xtrasl=alfa*x(ki)+beta
	  A1=cs*delta_x+xtrasl
	  B1=estremo
	 
	  alfa_j1=(B1-A1)/2.d0
	  beta_j1=(B1+A1)/2.d0

	  p2a = 0.d0
	  
	  IF ((B1-A1).gt.10.d-15) THEN

        DO kj=1,Ngauss

	        xinttrasl=(xint(kj)+1.d0)*0.5d0
	        s=fi1(iplog,iqlog,xinttrasl)
	        ds=dfi1(iplog,iqlog,xinttrasl)
	        serv=alfa_j1*(2.d0*s-1.d0)+beta_j1
            r2_1=(xtrasl-serv)**2

	        p2a = p2a + wint(kj)*ds*fiU(lj,serv,estremo,grado_q)*coeff_4*((cp**2*r2_1)**(-1))*cp*delta_x*dsqrt(dabs(-r2_1+cp**2*delta_x**2))
	    
	        IF (delta_kronecker(indice_i,indice_j).eq.1.d0) THEN
	    
	        p2a= p2a+wint(kj)*ds*fiU(lj,serv,estremo,grado_q)*coeff_5*(1.d0/(cp**2))*dlog(cp*delta_x+dsqrt(dabs(-r2_1+cp**2*delta_x**2)))
	        ENDIF
		  
        END DO

	  ENDIF

	  p2 = p2 +p2a*alfa_j1*w(ki)*fiU(li,xtrasl,estremo,grado_q)

    END DO
    Vudiag2_sub21= p2*alfa
    Vudiag_sub=Vudiag_sub+Vudiag2_sub21
    
  
!   primo integrale     
!   -----------------

	p3 = 0.d0
    
    IF (delta_kronecker(indice_i,indice_j).eq.1.d0) THEN
    
    DO ki=1,Ngauss
	    
	    p3a = 0.d0
	 
	    xtrasl=alfa*x(ki)+beta
	    A1=cs*delta_x+xtrasl
	    B1=estremo

	    alfa_j1=(B1-A1)/2.d0
	    beta_j1=(B1+A1)/2.d0
	  
	    if ((B1-A1).gt.10.d-15) then

		    cs1=(xtrasl-beta_j1)/alfa_j1

	        if (dabs(cs1).lt.1.25d0) then
                Write(*,*) 'ho usato la mopelog 2 cp'
				pause 
				
		        call mopelog(Ngaussi,ti,cs1,0.d0,ww1,1)
    	 
                DO kj=1,Ngaussi
			        p3a=p3a+vi(kj)*ww1(kj)*fiU(lj,alfa_j1*ti(kj)+beta_j1,estremo,grado_q)
                END DO

	            DO kj=1,Ngauss
	                p3a=p3a+w(kj)*fiU(lj,alfa_j1*x(kj)+beta_j1,estremo,grado_q)*dlog(alfa_j1)
                END DO

	        else
    	        
	            DO kj=1,Ngauss
	                p3a=p3a+w(kj)*fiU(lj,alfa_j1*x(kj)+beta_j1,estremo,grado_q)*(dlog(dabs(cs1-x(kj)))+dlog(alfa_j1))
                END DO
                
	        endif

	    endif

        p3=p3+p3a*alfa_j1*w(ki)*fiU(li,xtrasl,estremo,grado_q)
    END DO
    
    ENDIF  
  
    Vudiag2_sub22=-(1.d0/(cp**2))*coeff_5*p3*alfa
    Vudiag_sub=Vudiag_sub+Vudiag2_sub22
	
  	
  	!____________________________________intervallo 3__________________________
  	
  	
  	alfa=(sn-xn)/2.d0
	beta=(sn+xn)/2.d0

!   secondo integrale   
! ------------------
	iplog=2
    iqlog=1
	
	p2 = 0.d0
	 
	DO ki=1,Ngauss 
	
      xtrasl=alfa*x(ki)+beta
	  A1=0
	  B1=-cs*delta_x+xtrasl

	  alfa_j1=(B1-A1)/2.d0
	  beta_j1=(B1+A1)/2.d0

	  p2a = 0.d0
	  
	  IF ((B1-A1).gt.10.d-15) THEN

        DO kj=1,Ngauss

	        xinttrasl=(xint(kj)+1.d0)*0.5d0
	        s=fi1(iplog,iqlog,xinttrasl)
	        ds=dfi1(iplog,iqlog,xinttrasl)
	        serv=alfa_j1*(2.d0*s-1.d0)+beta_j1
            r2_1=(xtrasl-serv)**2

	        p2a = p2a + wint(kj)*ds*fiU(lj,serv,estremo,grado_q)*coeff_4*((cp**2*r2_1)**(-1))*cp*delta_x*dsqrt(dabs(-r2_1+cp**2*delta_x**2))
	    
	        IF (delta_kronecker(indice_i,indice_j).eq.1.d0) THEN
	    
	        p2a= p2a+wint(kj)*ds*fiU(lj,serv,estremo,grado_q)*coeff_5*(1.d0/(cp**2))*dlog(cp*delta_x+dsqrt(dabs(-r2_1+cp**2*delta_x**2)))
	        ENDIF
		  
        END DO

	  ENDIF

	  p2 = p2 +p2a*alfa_j1*w(ki)*fiU(li,xtrasl,estremo,grado_q)

    END DO
    Vudiag2_sub31= p2*alfa
    Vudiag_sub=Vudiag_sub+Vudiag2_sub31
   
  
!   primo integrale     
!   -----------------

	p3 = 0.d0
    
    IF (delta_kronecker(indice_i,indice_j).eq.1.d0) THEN
    
    DO ki=1,Ngauss
	    
	    p3a = 0.d0
	 
	    xtrasl=alfa*x(ki)+beta
	    A1=0
	    B1=-cs*delta_x+xtrasl

	    alfa_j1=(B1-A1)/2.d0
	    beta_j1=(B1+A1)/2.d0
	  
	    if ((B1-A1).gt.10.d-15) then

		    cs1=(xtrasl-beta_j1)/alfa_j1

	        if (dabs(cs1).lt.1.25d0) then
                Write(*,*) 'ho usato la mopelog 3 cp'
				pause
				
		        call mopelog(Ngaussi,ti,cs1,0.d0,ww1,1)
    	 
                DO kj=1,Ngaussi
			        p3a=p3a+vi(kj)*ww1(kj)*fiU(lj,alfa_j1*ti(kj)+beta_j1,estremo,grado_q)
                END DO

	            DO kj=1,Ngauss
	                p3a=p3a+w(kj)*fiU(lj,alfa_j1*x(kj)+beta_j1,estremo,grado_q)*dlog(alfa_j1)
                END DO

	        else
    	        
	            DO kj=1,Ngauss
	                p3a=p3a+w(kj)*fiU(lj,alfa_j1*x(kj)+beta_j1,estremo,grado_q)*(dlog(dabs(cs1-x(kj)))+dlog(alfa_j1))
                END DO
                
	        endif

	    endif

        p3=p3+p3a*alfa_j1*w(ki)*fiU(li,xtrasl,estremo,grado_q)
    END DO
    
    ENDIF  
  
    Vudiag2_sub32=-(1.d0/(cp**2))*coeff_5*p3*alfa
    Vudiag_sub=Vudiag_sub+Vudiag2_sub32
  endif
    
	!if((li.eq.1.and.lj.eq.2.or.lj.eq.1.and.li.eq.2))then!.and.delta_x.ge.0.1d0)then
	! write(*,*) 'integrale totale', Vudiag_sub
	! write(*,*) 'elemento', i
	! write(*,*) 'istante', delta_x
	! write(*,*) 'li', li
	! write(*,*) 'lj', lj
	! write(*,*) 'grado q', grado_q
	! write(*,*) 'integrali sul blu'
	! write(*,*) 'no log 1', Vudiag_sub11
	! write(*,*) 'log 1', Vudiag_sub12
	! write(*,*) 'no log 2', Vudiag_sub21
	! write(*,*) 'log 2', Vudiag_sub22
	! write(*,*) 'no log 2', Vudiag_sub31
	! write(*,*) 'log 2', Vudiag_sub32
	! write(*,*) 'integrali sul giallo'
    ! write(*,*) Vudiag2_sub11, Vudiag2_sub12
    ! write(*,*) Vudiag2_sub21, Vudiag2_sub22
    ! write(*,*) Vudiag2_sub31, Vudiag2_sub32
    ! write(*,*) Vudiag2_sub41, Vudiag2_sub42
	! write(*,*) 'integrale totale'
    ! write(*,*) Vudiag_sub
	!pause
	!endif
	


    RETURN
    END