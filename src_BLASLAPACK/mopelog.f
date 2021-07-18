	subroutine mopelog(n,t,as,bs,ww,ind)
c
c   modifica pesi per il nucleo log
c
	implicit real*8 (a-h,o-z)
	real*8 t(n),ww(n)

c       ind=1  ... elementi allineati (anche sovrapposti)
c       ind=2  ... elementi contigui non allineati
c
		if(ind.eq.1)then
	as0=1.d0+as
	as1=1.d0-as
	if(as1.eq.0.d0)then
	rmu0=as0*dlog(dabs(as0))-2.d0
	elseif(as0.eq.0.d0)then
	rmu0=as1*dlog(dabs(as1))-2.d0
	else
	rmu0=as0*dlog(dabs(as0))+as1*dlog(dabs(as1))-2.d0
	endif
c
	do 15 i=1,n
	ww(i)=rmu0
	do 10 k=1,n-1
	ww(i)=ww(i)+(2.d0*k+1)*rmu(k,as)*pleg(k,t(i))
10      continue

	ww(i)=ww(i)/2.d0
c
15      continue

		else

	as0=1.d0+as
	as1=1.d0-as
	if(as0.eq.0.d0.and.bs.eq.0.d0)then
	rmug0=2.d0*dlog(4.d0)-4.d0
	else
	rmug0=2.d0*bs*datan(as0/bs)+2.d0*bs*datan(as1/bs)+
     $        as0*dlog(bs**2+as0**2)+as1*dlog(bs**2+as1**2)-4.d0
	endif
c
	do 25 i=1,n
	ww(i)=rmug0
	do 20 k=1,n-1
	ww(i)=ww(i)+(2.d0*k+1)*rmug(k,as,bs)*pleg(k,t(i))
20      continue
c
	ww(i)=ww(i)/2.d0
25      continue
	      endif
	      return
	      end
c
c=======================================================================
c
       Real*8 function RMUG (kk,as,bs)

       Implicit double precision (a-h,o-z)

       a0=2.d0
       b0=2.d0/3.d0
       a1=15.d0/8.d0
       b1=3.d0/4.d0

       tm=4.d0/3.d0
       tm0=-as*tm
       tm1=b1*tm/a1
      tmm0=2.d0*bs*as*datan((as+1.d0)/bs)-2.d0*bs*as*datan((as-1.d0)/bs)
       tmm0=tmm0-0.5*(bs**2-as**2+1.d0)*dlog(bs**2+(as+1.d0)**2)
       tmm0=tmm0+0.5*(bs**2-as**2+1.d0)*dlog(bs**2+(as-1.d0)**2)-2.d0*as

       tmm1=2.d0*bs*(3.d0*as**2-(bs**2+1.d0))*(datan(as/bs+1.d0/bs)-
     ,datan(as/bs-1.d0/bs))+as*(as**2-(3.d0*bs**2+1.d0))*
     ,(dlog((as+1.d0)**2+bs**2)-dlog(bs**2+(as-1.d0)**2))-4.d0*
     ,(as**2-(bs**2+2.d0/3.d0))

       tmm2=a0*a1*(as*tmm1-(as**2+bs**2)*tmm0+tm0)-b1*tmm0

       q0i=bs*tmm0
       q0r=1.d0/a0*tmm1-as*tmm0
       rmug=q0i/(bs*kk)
	 if (kk.eq.1) RETURN
       q1i=bs*tmm1
       q1r=tmm2/a1+b1*tmm0/a1-as*tmm1
       rmug=q1i/(bs*kk)
	 if (kk.eq.2) RETURN
       q2i=(a1*as*q1i+a1*bs*q1r-b1*q0i)
       q2r=(a1*as*q1r-a1*bs*q1i-b1*q0r+a1*tm1)
       rmug=q2i/(bs*kk)
	 if (kk.eq.3) RETURN

       do 11, k=3,kk-1
       akk=((2*k+1.d0)*(k+1.d0))/(k*(k+2.d0))
       bkk=(k+1.d0)/(k+2.d0)

       qkki=akk*as*q2i+akk*bs*q2r-bkk*q1i
       qkkr=akk*as*q2r-akk*bs*q2i-bkk*q1r
       q1i=q2i
       q1r=q2r
       q2i=qkki
       q2r=qkkr
11     continue
       rmug=qkki/(bs*kk) 
       RETURN
       END
c
c=======================================================================
c
       Real*8 function pleg(kk,ss)
c
       implicit real*8 (a-h,o-z)
c
       p0=1.d0
       p1=ss
       pleg=p1
c
	 if (kk.eq.1) return
c
	   do 10 k=2,kk
	   coef=1.d0/k
	   p2=((2.d0*k-1)*ss*p1-(k-1)*p0)*coef
	   p0=p1
	   p1=p2
10         continue
c
       pleg=p2
c
       RETURN
       END
c
c==============================================================
c
	 Real*8 function rmu(kk,tt)
c
	 implicit real*8 (a-h,o-z)
c
	 if (dabs(tt).eq.1.d0) then
	  q110=-2.d0*tt
	 else
	  q110=(1.d0-tt**2)*dlog(dabs(1.d0-tt)/dabs(1.d0+tt))-2.d0*tt
	 endif
	 rmu=q110/2.d0
	 if (kk.eq.1) return
	 q111=2.d0*tt*q110+8.d0/3.d0
	 rmu=q111/4.d0
	 if(kk.eq.2)return
c
	     do 10  k=2,kk-1
		coef=(k+1.d0)/(k*(k+2.d0))
		q112=((2.d0*k+1.d0)*tt*q111-k*q110)*coef
	      q110=q111
	      q111=q112
10           continue
	 rmu=q112/(2.d0*kk)
c
	 RETURN
	 END
