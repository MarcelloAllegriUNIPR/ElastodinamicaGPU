      Double precision function Unota(j,x,y,t,indice_termine_noto)

      USE variable_2Dgeneral
  
      IMPLICIT NONE

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABILI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !Input
      INTEGER(kind=4),INTENT(IN):: j, indice_termine_noto
      REAL(kind=8),INTENT(IN):: x, y, t
  
      !Variabili locali
      INTEGER(kind=4),DIMENSION(2):: nodi
      REAl(kind=8),DIMENSION(2):: k
      
      REAL(kind=8):: elle, alfa, beta, coef, alfa_t, beta_t, F1_alfa, F2_beta, teta, tshift, omega, UnotaP, UnotaS, L, xbar, ybar
      
      !!!!!!!!!!!!!!!!!!!!!!!!!! CORPO della SUBROUTINE !!!!!!!!!!!!!!!!!!


    if (j.eq.1) then !termine noto che dà come soluzione 1
        nodi=list_elements(number_elements)%nodes
	    elle=list_nodes(nodi(2))%coordinates(1)
	    if((x-t).gt.0.d0) then
	        alfa=x-t
	    else 
	        alfa=0.d0
	    endif
	    if(elle.gt.(x+t)) then
	        beta=x+t
	    else
	        beta=elle
	    endif 
	    coef=1.d0/(8.d0*datan(1.d0))
	    if(t.eq.0.d0) then
            Unota= 0.d0
    	else
    	    alfa_t=x-alfa
    	    beta_t=beta-x
	        F1_alfa=-alfa_t*dlog(t+dsqrt(dabs((t-alfa_t)*(t+alfa_t))))+alfa_t*dlog(alfa_t)-t*datan(alfa_t/dsqrt(dabs((t-alfa_t)*(t+alfa_t))))
	        F2_beta=beta_t*dlog(t+dsqrt(dabs((t+beta_t)*(t-beta_t))))-beta_t*dlog(beta_t)+t*datan(beta_t/dsqrt(dabs((t+beta_t)*(t-beta_t))))
	        Unota=coef*(-F1_alfa+F2_beta)
     	endif
     	
    elseif (j.eq.2) then !plane linear wave Becache
        teta=4.d0*datan(1.d0)/2.d0
		tshift=t-dcos(teta)*x
		if (tshift.le.0.d0) then
			Unota=0.d0
		else
			Unota=tshift
		endif
    
    elseif (j.eq.3) then  !plane harmonic wave Becache
		teta=4.d0*datan(1.d0)*17.d0/36.d0   !4.d0*datan(1.d0)=pigreco
		omega=8.d0*(4.d0*datan(1.d0))
		tshift=t-dcos(teta)*x
		if (tshift.le.0.d0) then
			Unota=0.d0
		elseif (tshift.ge.0.d0.and.tshift.le.(4.d0*datan(1.d0)/omega)) then
			Unota=0.5d0*(1.d0-dcos(omega*tshift))
		else
			Unota=dsin(0.5d0*omega*tshift)
		endif
    
    elseif (j.eq.4) then   !impulso
		if ((t.le.(0.25d0+2.d0*deltat)).and.(t.ge.0.25d0)) then
			Unota=1.d0
		else
			Unota=0.d0
		endif
    
    elseif (j.eq.5) then   !tesi Brescia
		teta=4.d0*datan(1.d0)/2.d0   !4.d0*datan(1.d0)=pigreco
		omega=8.d0*(4.d0*datan(1.d0))
		tshift=t-dcos(teta)*x
		if (tshift.le.0.d0) then
			Unota=0.d0
		elseif (tshift.ge.0.d0.and.tshift.le.(4.d0*datan(1.d0)/omega)) then
			Unota=dsin(omega*tshift/2.d0)**2*x**4 !!!!!!!!!!!! togli x**4
		else
			Unota=1.d0*x**4 !!!!!!!!!!!!!!! x**4
		endif
    
    elseif (j.eq.6) then   !heaviside
		if (t.ge.0.0625d0) then
			Unota=1.d0
		else
			Unota=0.d0
		endif
    
    elseif (j.eq.7) then   !cappuccio
		if ((t.ge.0.25d0).and.(t.le.0.375d0)) then
			Unota=(t-0.25d0)/0.125d0
		elseif ((t.ge.0.375d0).and.(t.le.0.5d0)) then
			Unota=(-t+0.5d0)/0.125d0
		else
			Unota=0.d0
		endif
		
    elseif (j.eq.8) then   
		Unota=1.d0

	elseif(j.eq.9)then !!!!!!!!!!!!!AGGIUNTA ELASTODINAMICA!!!!!!!!!
	    if(indice_termine_noto.eq.1)then
	      if(t.lt.0.d0)then
	        Unota=0.d0
	      elseif(t.ge.0.d0.AND.t.le.(1.d0/8.d0))then
	        Unota=SIN(4.d0*(4.d0*datan(1.d0))*t)**2*x
	      else
	        Unota=x
	      endif
	    else
 	      if(t.lt.0.d0)then
	        Unota=0.d0
	      elseif(t.ge.0.d0.AND.t.le.(1.d0/8.d0))then
	        Unota=SIN(4.d0*(4.d0*datan(1.d0))*t)**2*dabs(x)
	      else
	        Unota=dabs(x)
	      endif	    
	    endif
	    
	elseif(j.eq.10)then !!!!!!!!!!!!!AGGIUNTA ELASTODINAMICA!!!!!!!!!
	    !if(indice_termine_noto.eq.1)then
	      !Unota=0.d0
	    !else
	      if(t.lt.0.d0)then
	        Unota=0.d0
	      elseif(t.ge.0.d0.AND.t.le.(1.d0/8.d0))then
	        Unota=SIN(4.d0*(4.d0*datan(1.d0))*t)**2*x
	      else
	        Unota=x
	      endif    
	    !endif
	elseif (j.eq.11) then   !prova simmetria
			Unota=3.d0
			
    elseif (j.eq.12) then 
         tshift=t
         if(tshift.le.0.5d0.or.tshift.gt.1.d0)then 
             Unota=0.d0
         elseif(tshift.gt.0.5d0.and.tshift.le.0.6d0)then
             Unota=tshift-0.5d0
         elseif(tshift.gt.0.6d0.and.tshift.le.0.7d0)then
             Unota=(-30.d0*tshift+19.d0)/10.d0
         elseif(tshift.gt.0.7d0.and.tshift.le.0.8d0)then
             Unota=(40.d0*tshift-30.d0)/10.d0
         elseif(tshift.gt.0.8d0.and.tshift.le.0.9d0)then
             Unota=(-30.d0*tshift+26.d0)/10.d0
         elseif(tshift.gt.0.9d0.and.tshift.le.1.d0)then
             Unota=tshift-1.d0
         endif
	
	 elseif (j.eq.13) then
         if(indice_termine_noto.eq.2)then	 
             tshift=t
             if(tshift.le.0.5d0.or.tshift.gt.1.d0)then 
                 Unota=0.d0
             elseif(tshift.gt.0.5d0.and.tshift.le.0.6d0)then
                 Unota=tshift-0.5d0
             elseif(tshift.gt.0.6d0.and.tshift.le.0.7d0)then
                 Unota=(-30.d0*tshift+19.d0)/10.d0
             elseif(tshift.gt.0.7d0.and.tshift.le.0.8d0)then
                 Unota=(40.d0*tshift-30.d0)/10.d0
             elseif(tshift.gt.0.8d0.and.tshift.le.0.9d0)then
                 Unota=(-30.d0*tshift+26.d0)/10.d0
             elseif(tshift.gt.0.9d0.and.tshift.le.1.d0)then
                 Unota=tshift-1.d0
             endif
		 else
		     Unota=0.d0
	     endif
	     
	 elseif(j.eq.14)then !!!!!!!!!!!!!AGGIUNTA ELASTODINAMICA!!!!!!!!!
	     if(t.lt.0.d0)then
		     Unota=0.d0
	     elseif(t.ge.0.d0.AND.t.le.(1.d0/8.d0))then
	         if(indice_termine_noto.eq.1)then
	             Unota=SIN(4.d0*(4.d0*datan(1.d0))*t)**2*x
	         else
	             Unota=SIN(4.d0*(4.d0*datan(1.d0))*t)**2*y
	         endif
	     else
	         if(indice_termine_noto.eq.1)then
	             Unota=x
	         else
	             Unota=y
	         endif
	     endif  
		 
	 elseif (j.eq.15) then	 
         tshift=t
         if(tshift.le.0.5d0.or.tshift.gt.1.d0)then 
             Unota=0.d0
         elseif(tshift.gt.0.5d0.and.tshift.le.0.6d0)then
             Unota=tshift-0.5d0
         elseif(tshift.gt.0.6d0.and.tshift.le.0.7d0)then
             Unota=(-30.d0*tshift+19.d0)/10.d0
         elseif(tshift.gt.0.7d0.and.tshift.le.0.8d0)then
             Unota=(40.d0*tshift-30.d0)/10.d0
         elseif(tshift.gt.0.8d0.and.tshift.le.0.9d0)then
             Unota=(-30.d0*tshift+26.d0)/10.d0
         elseif(tshift.gt.0.9d0.and.tshift.le.1.d0)then
             Unota=tshift-1.d0
         endif
		 if(indice_termine_noto.eq.1)then
             Unota=Unota*x
         else
             Unota=Unota*y		 
	     endif

	 elseif(j.eq.16)then !!!!!!!!!!!!!AGGIUNTA ELASTODINAMICA!!!!!!!!!
	     if(t.lt.0.d0.OR.t.gt.(1.d0/2.d0))then
		     Unota=0.d0
	     elseif(t.ge.0.d0.AND.t.le.(1.d0/8.d0))then
	         Unota=SIN(4.d0*(4.d0*datan(1.d0))*t)**2 
	     elseif(t.gt.(1.d0/8.d0).AND.t.lt.(3.d0/8.d0))then
	         Unota=1.d0
	     elseif(t.ge.(3.d0/8.d0).AND.t.le.(1.d0/2.d0))then
	         Unota=SIN(4.d0*(4.d0*datan(1.d0))*(1.d0/2.d0-t))**2
	     endif
	     
	     if(indice_termine_noto.eq.1)then
	             Unota=Unota*sign(1.d0,x)
	     else
	             !Unota=Unota*y
	             Unota=0.d0
	     endif 	 
	         
	 elseif(j.eq.17)then !!!!!!!!!!!!!AGGIUNTA ELASTODINAMICA!!!!!!!!!
	     if(t.lt.0.d0.OR.t.gt.0.48d0)then
		     Unota=0.d0
	     elseif(t.ge.0.d0.AND.t.le.0.12d0)then
	         Unota=SIN((25.d0/6.d0)*(4.d0*datan(1.d0))*t)**2 
	     elseif(t.gt.0.12d0.AND.t.lt.0.36d0)then
	         Unota=1.d0
	     elseif(t.ge.0.36d0.AND.t.le.0.48d0)then
	         Unota=SIN((25.d0/6.d0)*(4.d0*datan(1.d0))*(0.48d0-t))**2
	     endif
	     
	     if(indice_termine_noto.eq.1)then
	             !Unota=Unota*sign(1.d0,x)
	             Unota=Unota*x
	     else
	             Unota=Unota*y
	             !Unota=0.d0
	     endif 	   

	 elseif (j.eq.18) then	 
         tshift=t
         if(tshift.le.0.48d0.or.tshift.gt.0.88d0)then 
             Unota=0.d0
         elseif(tshift.gt.0.48d0.and.tshift.le.0.56d0)then
             Unota=tshift-0.48d0
         elseif(tshift.gt.0.56d0.and.tshift.le.0.64d0)then
             Unota=-3.d0*(tshift-0.56d0)+0.08d0
         elseif(tshift.gt.0.64d0.and.tshift.le.0.72d0)then
             Unota=4.d0*(tshift-0.64d0)-0.16d0;
         elseif(tshift.gt.0.72d0.and.tshift.le.0.8d0)then
             Unota=-3.d0*(tshift-0.72d0)+0.16d0;
         elseif(tshift.gt.0.8d0.and.tshift.le.0.88d0)then
             Unota=tshift-0.88d0
         endif
		 if(indice_termine_noto.eq.1)then
             Unota=Unota*x
         else
             Unota=Unota*y		 
	     endif   

	 elseif(j.eq.19)then !!!!!!!!!!!!!AGGIUNTA ELASTODINAMICA!!!!!!!!!
	     if(t.lt.0.d0)then
		     Unota=0.d0
	     elseif(t.ge.0.d0.AND.t.le.0.12d0)then
	         Unota=SIN((25.d0/6.d0)*(4.d0*datan(1.d0))*t)**2 
	     elseif(t.gt.0.12d0)then
	         Unota=1.d0
	     endif
	     
	     if(indice_termine_noto.eq.1)then
	             !Unota=Unota*sign(1.d0,x)
	             Unota=Unota*x
	     else
	             Unota=Unota*y
	             !Unota=0.d0
	     endif 	 

	 elseif(j.eq.20)then !!!!!!!!!!!!!AGGIUNTA ELASTODINAMICA!!!!!!!!!
         tshift=t
         if(tshift.le.0.5d0.or.tshift.gt.1.d0)then 
             Unota=0.d0
         elseif(tshift.gt.0.5d0.and.tshift.le.0.6d0)then
             Unota=tshift-0.5d0
         elseif(tshift.gt.0.6d0.and.tshift.le.0.7d0)then
             Unota=(-30.d0*tshift+19.d0)/10.d0
         elseif(tshift.gt.0.7d0.and.tshift.le.0.8d0)then
             Unota=(40.d0*tshift-30.d0)/10.d0
         elseif(tshift.gt.0.8d0.and.tshift.le.0.9d0)then
             Unota=(-30.d0*tshift+26.d0)/10.d0
         elseif(tshift.gt.0.9d0.and.tshift.le.1.d0)then
             Unota=tshift-1.d0
         endif
		 if(indice_termine_noto.eq.1)then
             Unota=Unota*y
         else
             Unota=Unota*(-x)		 
	     endif   
	    
	 elseif(j.eq.21)then !!!!!!!!!!!!!AGGIUNTA ELASTODINAMICA!!!!!!!!!
         tshift=t
         if(tshift.le.0.5d0.or.tshift.gt.1.d0)then 
             Unota=0.d0
         elseif(tshift.gt.0.5d0.and.tshift.le.0.6d0)then
             Unota=tshift-0.5d0
         elseif(tshift.gt.0.6d0.and.tshift.le.0.7d0)then
             Unota=(-30.d0*tshift+19.d0)/10.d0
         elseif(tshift.gt.0.7d0.and.tshift.le.0.8d0)then
             Unota=(40.d0*tshift-30.d0)/10.d0
         elseif(tshift.gt.0.8d0.and.tshift.le.0.9d0)then
             Unota=(-30.d0*tshift+26.d0)/10.d0
         elseif(tshift.gt.0.9d0.and.tshift.le.1.d0)then
             Unota=tshift-1.d0
         endif
         teta=datan(-1.d0)
		 if(indice_termine_noto.eq.1)then
             Unota=Unota*(x*dcos(teta)-y*dsin(teta))
         else
             Unota=Unota*(x*dsin(teta)+y*dcos(teta))	 
	     endif   

	 elseif(j.eq.22)then !!!!!!!!!!!!!AGGIUNTA ELASTODINAMICA!!!!!!!!!
         tshift=t-1
         if(tshift.le.0.5d0.or.tshift.gt.1.d0)then 
             Unota=0.d0
         elseif(tshift.gt.0.5d0.and.tshift.le.0.6d0)then
             Unota=tshift-0.5d0
         elseif(tshift.gt.0.6d0.and.tshift.le.0.7d0)then
             Unota=(-30.d0*tshift+19.d0)/10.d0
         elseif(tshift.gt.0.7d0.and.tshift.le.0.8d0)then
             Unota=(40.d0*tshift-30.d0)/10.d0
         elseif(tshift.gt.0.8d0.and.tshift.le.0.9d0)then
             Unota=(-30.d0*tshift+26.d0)/10.d0
         elseif(tshift.gt.0.9d0.and.tshift.le.1.d0)then
             Unota=tshift-1.d0
         endif
		 if(indice_termine_noto.eq.1)then
             Unota=Unota*x
         else
             Unota=Unota*y
	     endif  

	 elseif(j.eq.23)then !!!!!!!!!!!!!AGGIUNTA ELASTODINAMICA!!!!!!!!!
         !tshift=velC_P*(t-0.5d0)-x
         tshift=2.5d0*(t-0.5d0)-x
         if(tshift.le.0.5d0.or.tshift.gt.1.d0)then 
             UnotaP=0.d0
         elseif(tshift.gt.0.5d0.and.tshift.le.0.6d0)then
             UnotaP=tshift-0.5d0			 
         elseif(tshift.gt.0.6d0.and.tshift.le.0.7d0)then
             UnotaP=(-30.d0*tshift+19.d0)/10.d0
         elseif(tshift.gt.0.7d0.and.tshift.le.0.8d0)then
             UnotaP=(40.d0*tshift-30.d0)/10.d0
         elseif(tshift.gt.0.8d0.and.tshift.le.0.9d0)then
             UnotaP=(-30.d0*tshift+26.d0)/10.d0
         elseif(tshift.gt.0.9d0.and.tshift.le.1.d0)then
             UnotaP=tshift-1.d0
         endif
		 tshift=velC_S*(t-0.5d0)-x
         if(tshift.le.0.5d0.or.tshift.gt.1.d0)then 
             UnotaS=0.d0
         elseif(tshift.gt.0.5d0.and.tshift.le.0.6d0)then
             UnotaS=tshift-0.5d0			 
         elseif(tshift.gt.0.6d0.and.tshift.le.0.7d0)then
             UnotaS=(-30.d0*tshift+19.d0)/10.d0
         elseif(tshift.gt.0.7d0.and.tshift.le.0.8d0)then
             UnotaS=(40.d0*tshift-30.d0)/10.d0
         elseif(tshift.gt.0.8d0.and.tshift.le.0.9d0)then
             UnotaS=(-30.d0*tshift+26.d0)/10.d0
         elseif(tshift.gt.0.9d0.and.tshift.le.1.d0)then
             UnotaS=tshift-1.d0
         endif
		 if(indice_termine_noto.eq.1)then
             Unota=-UnotaP
		 else
             Unota=0.d0 
	     endif   		 
 
 	 elseif(j.eq.24)then !!!!!!!!!!!!!AGGIUNTA ELASTODINAMICA!!!!!!!!!
	     if(t.le.0.d0)then
	         Unota=0.d0
		 elseif(t.gt.0.d0.AND.t.le.1.d0/velC_P)then
		     Unota=(1.d0/velC_P**2)*(1.d0/(velC_P**2*(t+1.d0/velC_P)**2-1)**(3.d0/2.d0))
		 else
		     Unota=0.d0
	     endif
		 if(indice_termine_noto.eq.1)then
             Unota=Unota*x
		 else
             Unota=Unota*y
	     endif 

 	 elseif(j.eq.25)then !!!!!!!!!!!!!AGGIUNTA ELASTODINAMICA!!!!!!!!!
	     if(t.le.0.d0)then
	         Unota=0.d0
		 elseif(t.gt.0.d0.AND.t.le.1.d0/velC_S)then
		     Unota=(1.d0/30.d0)*(1.d0/(velC_S**2*(t+1.d0/velC_S)**2-1)**(3.d0/2.d0))
		 else
		     Unota=0.d0
	     endif
		 if(indice_termine_noto.eq.1)then
             Unota=Unota*y
         else
             Unota=Unota*(-x)		 
	     endif    	

 	 elseif(j.eq.26)then !!!!!!!!!!!!!AGGIUNTA ELASTODINAMICA!!!!!!!!!
	     if(t.le.0.d0)then
	         UnotaP=0.d0
		 elseif(t.gt.0.d0.AND.t.le.1.d0/velC_P)then
		     UnotaP=(velC_P/(2.d0*4.d0*datan(1.d0)*10.d0))*(1.d0/(velC_P**2*(t+1.d0/velC_P)**2-1)**(3.d0/2.d0))
		 else
		     UnotaP=0.d0
	     endif
		 if(t.le.0.d0)then
	         UnotaS=0.d0
		 elseif(t.gt.0.d0.AND.t.le.1.d0/velC_S)then
		     UnotaS=(velC_S/(2.d0*4.d0*datan(1.d0)*10.d0))*(1.d0/(velC_S**2*(t+1.d0/velC_S)**2-1)**(3.d0/2.d0))
		 else
		     UnotaS=0.d0
	     endif 
		 if(indice_termine_noto.eq.1)then
             Unota=UnotaP*x+UnotaS*y
		 else
             Unota=UnotaP*y+UnotaS*(-x)
	     endif 	 		 

	 elseif(j.eq.27)then !!!!!!!!!!!!!AGGIUNTA ELASTODINAMICA!!!!!!!!!
		 tshift=velC_S*(t-0.75d0)-x
         if(tshift.le.0.5d0.or.tshift.gt.1.d0)then 
             UnotaS=0.d0
         elseif(tshift.gt.0.5d0.and.tshift.le.0.6d0)then
             UnotaS=tshift-0.5d0			 
         elseif(tshift.gt.0.6d0.and.tshift.le.0.7d0)then
             UnotaS=(-30.d0*tshift+19.d0)/10.d0
         elseif(tshift.gt.0.7d0.and.tshift.le.0.8d0)then
             UnotaS=(40.d0*tshift-30.d0)/10.d0
         elseif(tshift.gt.0.8d0.and.tshift.le.0.9d0)then
             UnotaS=(-30.d0*tshift+26.d0)/10.d0
         elseif(tshift.gt.0.9d0.and.tshift.le.1.d0)then
             UnotaS=tshift-1.d0
         endif
		 if(indice_termine_noto.eq.1)then
             Unota=0.d0 
		 else
             Unota=UnotaS
	     endif   		 

 	 elseif(j.eq.28)then !!!!!!!!!!!!!AGGIUNTA ELASTODINAMICA!!!!!!!!!
	     !if(t.le.0.d0)then
	      !   UnotaP=0.d0
		 if(t.gt.0.d0)then!.AND.t.le.2.d0/velC_P)then
		     UnotaP=(velC_P*50.d0/(2.d0*4.d0*datan(1.d0)))*(1.d0/(velC_P**2*(t+2.d0/velC_P)**2-1)**(3.d0/2.d0))
		 else
		     UnotaP=0.d0
	     endif
		 !if(t.le.0.d0)then
	      !   UnotaS=0.d0
		 if(t.gt.0.d0)then!.AND.t.le.2.d0/velC_S)then
		     UnotaS=(velC_S*50.d0/(2.d0*4.d0*datan(1.d0)))*(1.d0/(velC_S**2*(t+2.d0/velC_S)**2-1)**(3.d0/2.d0))
		 else
		     UnotaS=0.d0
	     endif 
		 if(indice_termine_noto.eq.1)then
             Unota=UnotaP*x+UnotaS*y
		 else
             Unota=UnotaP*y+UnotaS*(-x)
	     endif	 
		 
		     
 	 elseif(j.eq.29)then !!!!!!!!!!!!!AGGIUNTA ELASTODINAMICA!!!!!!!!!
 	     !k(1)=dcos(4.d0*datan(1.d0)/3.d0)
 	     !k(2)=dsin(4.d0*datan(1.d0)/3.d0)
		 !k(1)=-dcos(datan(1.d0))
		 !k(2)=dsin(datan(1.d0))
		 k(1)=0.d0
		 k(2)=1.d0
		 tshift=velC_P*(t-0.5d0)!-k(1)*x-k(2)*y
         if(tshift.le.0.5d0.or.tshift.gt.1.d0)then 
             Unota=0.d0
         elseif(tshift.gt.0.5d0.and.tshift.le.0.6d0)then
             Unota=tshift-0.5d0			 
         elseif(tshift.gt.0.6d0.and.tshift.le.0.7d0)then
             Unota=(-30.d0*tshift+19.d0)/10.d0
         elseif(tshift.gt.0.7d0.and.tshift.le.0.8d0)then
             Unota=(40.d0*tshift-30.d0)/10.d0
         elseif(tshift.gt.0.8d0.and.tshift.le.0.9d0)then
             Unota=(-30.d0*tshift+26.d0)/10.d0
         elseif(tshift.gt.0.9d0.and.tshift.le.1.d0)then
             Unota=tshift-1.d0
         endif
		 if(indice_termine_noto.eq.1)then
             Unota=-k(1)*Unota
		 else
             Unota=-k(2)*Unota
	     endif 
	     
	elseif(j.eq.30)then !!!!!!!!!!!!!AGGIUNTA ELASTODINAMICA!!!!!!!!!
	     if(t.le.0.d0)then
	         UnotaP=0.d0
		 elseif(t.gt.0.d0)then !.AND.t.le.1.d0/velC_P)then
		     UnotaP=(velC_P/(2.d0*4.d0*datan(1.d0)))*(1.d0/(velC_P**2*(t+1.d0/velC_P)**2-1)**(1.d0/2.d0))
		 !else
		     !UnotaP=0.d0
	     endif
		 if(t.le.0.d0)then
	         UnotaS=0.d0
		 elseif(t.gt.0.d0)then !.AND.t.le.1.d0/velC_S)then
		     UnotaS=(velC_S/(2.d0*4.d0*datan(1.d0)))*(1.d0/(velC_S**2*(t+1.d0/velC_S)**2-1)**(1.d0/2.d0))
		 !else
		     !UnotaS=0.d0
	     endif 
		 if(indice_termine_noto.eq.1)then
             Unota=UnotaP*x+UnotaS*y
		 else
             Unota=UnotaP*y+UnotaS*(-x)
	     endif 
	     
	 elseif(j.eq.31)then !!!!!!!!!!!!!AGGIUNTA ELASTODINAMICA!!!!!!!!!
	    if(indice_termine_noto.eq.2)then
	      if(t.lt.0.d0)then
	        Unota=0.d0
	      elseif(t.ge.0.d0.AND.t.le.(1.d0/8.d0))then
	        Unota=SIN(4.d0*(4.d0*datan(1.d0))*t)**2*dabs(x)**10!*((0.5d0**2-x**2)**(0.25)+0.5d0)
	      else
	        Unota=1*dabs(x)**10!*((0.5d0**2-x**2)**(0.25)+0.5d0)
	      endif
	    else
          !if(t.lt.0.d0)then
	        Unota=0.d0
	      !elseif(t.ge.0.d0.AND.t.le.(1.d0/8.d0))then
	        !Unota=SIN(4.d0*(4.d0*datan(1.d0))*t)**2*x!**8 !*((3.d0/4.d0-y**2)**(0.25)+0.5d0)
	      !else
	        !Unota=1*x !**8 !*((3.d0/4.d0-y**2)**(0.25)+0.5d0)
	      !endif
	    endif	 
	    
	   elseif(j.eq.32)then !!!!!!!!!!!!!AGGIUNTA ELASTODINAMICA!!!!!!!!!
	    if(indice_termine_noto.eq.1)then
	      Unota=0.d0
	    else
	      if((t-y).lt.0.d0)then
	        Unota=0.d0
	      elseif((t-y).ge.0.d0.AND.t.le.(1.d0/8.d0))then
	        Unota=SIN(4.d0*(4.d0*datan(1.d0))*(t-y))**2
	      else
	        Unota=1.d0
	      endif    
	    endif		         
	   
	  elseif(j.eq.33)then !!!!!!!!!!!!!AGGIUNTA ELASTODINAMICA!!!!!!!!!
	    if(indice_termine_noto.eq.1)then
	      Unota=0.d0
	    else
		  beta=4.d0*datan(1.d0)/3.d0 !pi/3 !angolo top
		  !beta=4.d0*datan(1.d0)/12.d0 !pi/3 !angolo top
		  !beta=3.d0*datan(1.d0) !pi/3 !angolo top
		  alfa=(4.d0*datan(1.d0)-beta)/2.d0 !pi/3 !angolo base
		  L=dcos(alfa) !semilunghezza base
	      if((t+y-L*dtan(alfa)).lt.0.d0)then
	        Unota=0.d0
	      elseif((t+y-L*dtan(alfa)).ge.0.d0.AND.t.le.(1.d0/8.d0))then
	        Unota=SIN(4.d0*(4.d0*datan(1.d0))*(t+y-L*dtan(alfa)))**2
	      else
	        Unota=1.d0
	      endif    
	    endif	
	    
	  elseif(j.eq.34)then !!!!!!!!!!!!!AGGIUNTA ELASTODINAMICA!!!!!!!!!
	    if(indice_termine_noto.eq.1)then
	      Unota=0.d0
	    else
	      if(t.lt.0.d0)then
	        Unota=0.d0
	      elseif(t.ge.0.d0.AND.t.le.(1.d0/8.d0))then
	        Unota=SIN(4.d0*(4.d0*datan(1.d0))*t)**2*dabs(y-dsqrt(3.d0)/4.d0)**4.d0 !**(1.d0/4.d0)
	      else
	        Unota=dabs(y-dsqrt(3.d0)/4.d0)**4.d0 !!**(1.d0/4.d0)
	      endif    
	    endif	
		
	  elseif(j.eq.35)then !!!!!!!!!!!!!AGGIUNTA ELASTODINAMICA!!!!!!!!!
	     !write(*,*) 'sono quiiii!!!!!!!!!'
		 xbar=0.d0
		 !ybar=0.127561144121697d0
		 ybar=1.d0/(2.d0*dsqrt(3.d0))
		 if(t.lt.0.d0)then
	        Unota=0.d0
	     elseif(t.ge.0.d0.AND.t.le.(1.d0/8.d0))then
	        Unota=dsin(4.d0*(4.d0*datan(1.d0))*t)**2
	     else
	        Unota=1.d0
	     endif
		 ! if(y.eq.0.d0)then
		     ! k(1)=0.d0
			 ! k(2)=-1.d0
		 ! else
		     ! if(x.le.xbar)then
			     ! k(1)=dcos(5.d0*4.d0*datan(1.d0)/8.d0)
			     ! k(2)=dsin(5.d0*4.d0*datan(1.d0)/8.d0)
			 ! else
			     ! k(1)=dcos(3.d0*4.d0*datan(1.d0)/8.d0)
			     ! k(2)=dsin(3.d0*4.d0*datan(1.d0)/8.d0)
			 ! endif
		 ! endif
		 if(y.eq.(ybar-1.d0/(2.d0*dsqrt(3.d0))))then
		     k(1)=0
			 k(2)=-1.d0
		 else
		     if(x.le.xbar)then
			     k(1)=dcos(5.d0*4.d0*datan(1.d0)/6.d0)
			     k(2)=dsin(5.d0*4.d0*datan(1.d0)/6.d0)
			 else
			     k(1)=dcos(4.d0*datan(1.d0)/6.d0)
			     k(2)=dsin(4.d0*datan(1.d0)/6.d0)
			 endif
		 endif
		 if(indice_termine_noto.eq.1)then
             Unota=k(1)*Unota
		 else
             Unota=k(2)*Unota
	     endif
		 
	   elseif(j.eq.36)then !!!!!!!!!!!!!AGGIUNTA ELASTODINAMICA!!!!!!!!!
	      if(t.lt.0.d0)then
	        Unota=0.d0
	      elseif(t.ge.0.d0.AND.t.le.(1.d0/8.d0))then
	        Unota=SIN(4.d0*(4.d0*datan(1.d0))*t)**2 !*dabs(dabs(x)-dcos(datan(1.d0)/2.d0)) !**4
	      else
	        Unota=1 !*dabs(dabs(x)-dcos(datan(1.d0)/2.d0)) !**4
	      endif 
		     k(1)=-dcos(5.d0*4.d0*datan(1.d0)/8.d0)
			 k(2)=-dsin(5.d0*4.d0*datan(1.d0)/8.d0)
			 !k(1)=0.d0
			 !k(2)=1.d0
			 !k(1)=-dcos(3.d0*datan(1.d0))
			 !k(2)=-dsin(3.d0*datan(1.d0))
          if(indice_termine_noto.eq.1)then
             Unota=k(1)*Unota
		  else
             Unota=k(2)*Unota
	      endif		  
	 
		 
      else
	 
	    write (*,*) 'Attenzione!Funzione non definita entro il set.'
        write(*,*)'sono in Unota',j
        pause
        stop
     endif   

    RETURN 
    END
	
	!xbar=0.d0
	!ybar=1.d0/(2.d0*dsqrt(3.d0))
	!xbar=0.d0
	!ybar=0.d0
	 
		 
		 ! if(y.eq.(ybar-1.d0/(2.d0*dsqrt(3.d0))))then
		     ! k(1)=0
			 ! k(2)=-1.d0
		 ! else
		     ! if(x.le.xbar)then
			     ! k(1)=dcos(5.d0*4.d0*datan(1.d0)/6.d0)
			     ! k(2)=dsin(5.d0*4.d0*datan(1.d0)/6.d0)
			 ! else
			     ! k(1)=dcos(4.d0*datan(1.d0)/6.d0)
			     ! k(2)=dsin(4.d0*datan(1.d0)/6.d0)
			 ! endif
		 ! endif
		 ! if(dabs(x-(xbar+1.d0/(2.d0*dsqrt(3.d0)))).le.1.d-12)then !triangolo ruotato di 45gradi
		     ! k(1)=1.d0
			 ! k(2)=0.d0
		 ! else
		     ! if(y.ge.xbar)then
			     ! k(1)=dcos(2.d0*4.d0*datan(1.d0)/3.d0)
			     ! k(2)=dsin(2.d0*4.d0*datan(1.d0)/3.d0)
			 ! else
			     ! k(1)=dcos(-2.d0*4.d0*datan(1.d0)/3.d0)
			     ! k(2)=dsin(-2.d0*4.d0*datan(1.d0)/3.d0)
			 ! endif
		 ! endif