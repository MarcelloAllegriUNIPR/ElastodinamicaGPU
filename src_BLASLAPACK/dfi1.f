	real*8 function dfi1(ip,iq,t)

	implicit double precision (a-h,o-z)

c       coefficiente
	cn=ip
	cd=1.d0
	do 15 i=1,iq-1
	cd=cd*i
	cn=cn*(ip+i)
15      continue
	coef=cn/cd

	dfi1=coef*t**(ip-1)*(1.d0-t)**(iq-1)

	return
	end
