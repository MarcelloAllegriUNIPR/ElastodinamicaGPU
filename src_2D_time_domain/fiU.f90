      Double precision function fiU(l,x,lung,grado_q)   
!      
!     FUNZIONE DI FORMA L-ESIMA SULL'ELEMENTO i, SU UN SET DI LU
! 
  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABILI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  INTEGER(kind=4),INTENT(IN):: l, grado_q
  REAL(kind=8),INTENT(IN):: x, lung
  
  !Variabili locali
  REAL(kind=8),EXTERNAL::plagran
      
      if ((l.lt.1).or.(l.gt.(grado_q+1))) then
       write (*,*) 'Errore di indice nella funzione di forma di gammau'
      else
       fiU=plagran(grado_q,l-1,x/(lung/2.d0)-1.d0)
      endif

  RETURN
  END
