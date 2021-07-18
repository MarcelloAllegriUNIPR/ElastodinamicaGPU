      Double precision function plagran(igrado,j,x)   
                                           !Funzione che valuta 
                                           !il polinomio j-esimo
                                           !interpolatore di Lagrange 
                                           !di grado igrado, su una 
                                           !decomposizione uniforme 
                                           !dell'intervallo, nel punto x
      implicit double precision (a-h,o-z)

      pj=1.d0
      if (igrado.gt.0) then
!        do 10,k=0,igrado
         DO k=0,igrado
            if (k.ne.j) then
              pnum=igrado*(x+1)-2*k
              den=2*(j-k)
              pj=pj*pnum/den
              endif
!10       continue 
         END DO
      endif
      plagran=pj
      RETURN
      END
