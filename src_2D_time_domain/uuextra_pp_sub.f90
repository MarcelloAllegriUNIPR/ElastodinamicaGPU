Double precision function uuextra_pp_sub(j,lj,x,y,tempo1,velC,velP,tempo2,indice_termine_noto,indice_j)

!
!     INTEGRALE CON FUNZIONE DI FORMA lj SULL'ELEMENTO j
!
  
  USE variable_2Dgeneral
  
  IMPLICIT NONE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABILI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  INTEGER(kind=4),INTENT(IN):: j, lj, indice_termine_noto, indice_j
  REAL(kind=8),INTENT(IN)::velC, velP, x, y, tempo1, tempo2
  
  !Variabili locali
  INTEGER(kind=4):: AllocateStatus, ind_gauss, Ngauss, kj, k, l, jj
  INTEGER(kind=4),DIMENSION(2):: nodi_j
  REAL(kind=8),DIMENSION(:),ALLOCATABLE:: xx, ww
  REAL(kind=8),DIMENSION(2):: coor_j1, coor_j2, r_uuextra
  REAL(kind=8),DIMENSION(4):: p, DatoP, DatoS
  REAL(kind=8):: deltax, deltay, beta, gammaS, gammaP, dxS, dxP, xpos, ypos, r2, & 
                 len_j, csi1, csi2, cpi1, cpi2, app, dx, delta, dintP, dintS, &
				 xxP, wwP, xxS, wwS, minimo, massimo, rad
  REAL(kind=8),EXTERNAL:: fiU
  INTEGER(kind=4), DIMENSION(2,2) :: delta_kronecker

  !!!!!!!!!!!!!!!!!!!!!!!!!! CORPO della SUBROUTINE !!!!!!!!!!!!!!!!!!
	delta_kronecker(1,1)=1
    delta_kronecker(1,2)=0
    delta_kronecker(2,1)=0
    delta_kronecker(2,2)=1

	
	uuextra_pp_sub=0.d0
	
	if ((tempo1-tempo2).le.0.d0) return
	
	ind_gauss=6
    Ngauss=2**(ind_gauss-1) !32 nodi di Gauss
	
    ALLOCATE(xx(Ngauss),STAT=AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    ALLOCATE(ww(Ngauss),STAT=AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
	
    xx=gauss(ind_gauss)%nodiquad
    ww=gauss(ind_gauss)%pesiquad
    
    !Estremi sottointervallo j 
    nodi_j=list_elements(j)%nodes

    !Coordinate degli estremi 
    coor_j1=list_nodes(nodi_j(1))%coordinates(1:2)
    coor_j2=list_nodes(nodi_j(2))%coordinates(1:2)
    
    deltax=coor_j2(1)-coor_j1(1)
    deltay=coor_j2(2)-coor_j1(2)
	len_j=(list_elements(j)%length)
	
	dxS=velC*(tempo1-tempo2)
	dxP=velP*(tempo1-tempo2)
	dx=tempo1-tempo2
	
	beta=2.d0*(deltax*(coor_j1(1)-x)+deltay*(coor_j1(2)-y))/len_j
	gammaS=(x-coor_j1(1))**2+(y-coor_j1(2))**2-dxS**2
	gammaP=(x-coor_j1(1))**2+(y-coor_j1(2))**2-dxP**2
	
	if ((beta**2-4.d0*gammaP).le.(10.d0**(-13))) return
	
	cpi1=(-beta-dsqrt(beta**2-4.d0*gammaP))/2.d0
	cpi2=(-beta+dsqrt(beta**2-4.d0*gammaP))/2.d0
    
    deltax=deltax/len_j
    deltay=deltay/len_j
	
	jj=indice_j
    
	!integraione nucleo con attiva la velocità velP
    dintP=0.d0
	minimo=dmin1(len_j,cpi2)
	massimo=dmax1(0.d0,cpi1)
	delta=minimo-massimo
	if (delta.gt.10.d0**(-11))then
         DO kj=1,ngauss
	         xxP=delta*0.5d0*(xx(kj)+1.d0)+massimo
	         wwP=ww(kj)/2.d0
	         xpos=coor_j1(1)+xxP*deltax
	         ypos=coor_j1(2)+xxP*deltay
			 r_uuextra(1)=x-xpos
			 r_uuextra(2)=y-ypos
	         r2=r_uuextra(1)**2+r_uuextra(2)**2   
			 rad=dmax1(dxP**2-r2,0.d0)
			 !rad=dxP**2-r2
			 
			 dintP=dintP+wwP*fiU(lj,xxP,len_j,grado_q)*&
			 !(dx*dsqrt(dxP**2-r2)/velP)
			 ((r_uuextra(jj)*r_uuextra(indice_termine_noto)/r2**2-delta_kronecker(jj,indice_termine_noto)/(2*r2))*dx*dsqrt(rad)/velP+&
			 (delta_kronecker(jj,indice_termine_noto)/2.d0)*(1/velP**2)*(dlog(dxP+dsqrt(rad))-dlog(dsqrt(r2))))
			 ! if ((dxP**2-r2).lt.0.d0)then
			     ! write(*,*) j, lj, x, y
				 ! write(*,*) velC, velP, tempo1, tempo2
				 ! write(*,*) indice_termine_noto, jj
				 ! write(*,*) dxP**2-r2 
				 ! pause
			 ! endif
         END DO
    endif
    dintP=dintP*delta
	
	!integraione nucleo con attive entrambe solo la velC	
    dintS=0.d0
	if ((beta**2-4.d0*gammaS).ge.(10.d0**(-13)))then
	     csi1=(-beta-dsqrt(beta**2-4.d0*gammaS))/2.d0
	     csi2=(-beta+dsqrt(beta**2-4.d0*gammaS))/2.d0

       	 minimo=dmin1(len_j,csi2)
	     massimo=dmax1(0.d0,csi1)
         delta=minimo-massimo
	     if (delta.gt.10.d0**(-11))then 
             DO kj=1,ngauss
	             xxS=delta*0.5d0*(xx(kj)+1.d0)+massimo
	             wwS=ww(kj)/2.d0
	             xpos=coor_j1(1)+xxS*deltax
	             ypos=coor_j1(2)+xxS*deltay
			     r_uuextra(1)=x-xpos
			     r_uuextra(2)=y-ypos
	             r2=r_uuextra(1)**2+r_uuextra(2)**2   
			     rad=dmax1(dxS**2-r2,0.d0)
				 !rad=dxs**2-r2
				 
			     dintS=dintS+wwS*fiU(lj,xxS,len_j,grado_q)*&
				 !(r_uuextra(jj)*r_uuextra(indice_termine_noto)/r2**2)
			     (-(r_uuextra(jj)*r_uuextra(indice_termine_noto)/r2**2-delta_kronecker(jj,indice_termine_noto)/(2*r2))*dx*dsqrt(rad)/velC+&
			     (delta_kronecker(jj,indice_termine_noto)/2.d0)*(1/velC**2)*(dlog(dxS+dsqrt(rad))-dlog(dsqrt(r2))))
		     END DO
         endif
    endif
    dintS=dintS*delta	
	
    uuextra_pp_sub=(1.d0/(8.d0*datan(1.d0)*rho))*(dintP+dintS)
    
RETURN
END
