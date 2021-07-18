	real*8 function fi1(ip,iq,t)

	implicit double precision (a-h,o-z)
	parameter (m=64)

	dimension bser(m),endpts(2),fn(m),fw(m)

	Ngauss=(ip+iq)/2

	call gaussq(1,Ngauss,0.d0,0.d0,0,endpts,bser,fn,fw)

c       traslo i nodi

	do 10 i=1,Ngauss
	fn(i)=0.5d0*t*(fn(i)+1.d0)
	fw(i)=0.5d0*t*fw(i)
10      continue

c       coefficiente
	cn=ip
	cd=1.d0
	do 15 i=1,iq-1
	cd=cd*i
	cn=cn*(ip+i)
15      continue
	coef=cn/cd

	fi1=0.d0
	do 25i=1,Ngauss
	fi1=fi1+fw(i)*fn(i)**(ip-1)*(1.d0-fn(i))**(iq-1)
25      continue
	fi1=fi1*coef

	return
	end
