  logical function ULeftAss(i)
                                      !Funzione che verifica se l'elemento i
                                      !e` "assemblato con un elemento sulla sua sinistra".
  USE variable_2Dgeneral
  
  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABILI !!!!!!!!!!!!!!!!!!!!!!!!!
  
  !Input
  INTEGER(kind=4),INTENT(IN):: i

  !Variabili locali
  INTEGER(kind=4):: j

  !!!!!!!!!!!!!!!!!!!!!!!!!! CORPO della SUBROUTINE !!!!!!!!!!!!!!!!!!
  
  DO j=1,number_elements
    IF ((list_elements(i)%nodes(1).EQ.list_elements(j)%nodes(2)).and.&
    (list_elements(i)%normal(1).EQ.list_elements(j)%normal(1)).and.&
    (list_elements(i)%normal(2).EQ.list_elements(j)%normal(2)).and.&
    (list_elements(i)%normal(3).EQ.list_elements(j)%normal(3))) THEN
        ULeftAss=.true.
        exit
    ELSE
        ULeftAss=.false.
    END IF
  END DO
  
  RETURN
  END