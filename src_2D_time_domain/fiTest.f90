      Double precision function fiTest(i,l,x)
!      
!     FUNZIONE DI FORMA L-ESIMA SULL'ELEMENTO i, SU UN SET DI LQ
!
      USE variable_2Dgeneral
  
      IMPLICIT NONE

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABILI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !Input
      INTEGER(kind=4),INTENT(IN):: i, l
      REAL(kind=8),INTENT(IN):: x
  
      !Variabili locali
      REAL(kind=8):: len_i
      integer :: ltest
      
      REAL(kind=8),EXTERNAL::plagran
      
      !!!!!!!!!!!!!!!!!!!!!!!!!! CORPO della SUBROUTINE !!!!!!!!!!!!!!!!!!

      ltest=grado_q+1
      len_i=(list_elements(i)%length)/2.d0
	  if ((l.lt.1).or.(l.gt.ltest)) then
        write (*,*) 'Errore di indice nella funzione di forma di gammaq'
      else
        fiTest=plagran(ltest-1,l-1,x/len_i-1.d0)
      endif

      RETURN
      END