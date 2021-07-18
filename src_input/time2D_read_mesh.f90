SUBROUTINE time2D_read_mesh(filemesh)

  USE variable_2Dgeneral

  IMPLICIT NONE

  !Input 
  CHARACTER(len=200),INTENT(IN):: filemesh

  !Variabili locali
  CHARACTER(len=200):: rec

  !!!!!!!!!!!!!!!!!!!!!!!!!! CORPO della SUBROUTINE !!!!!!!!!!!!!!!!!!

  !Apertura del file di input
  OPEN(10,FILE=filemesh,STATUS='OLD')
  !Ci posizioniamo all'inizio del file di input
  REWIND(10)! fa partire l'istruzione READ dalla prima riga del file

  !Lettura dell'intestazione del file mesh
  READ(10,*) rec
  DO WHILE(rec(1:8).NE.'Vertices')
     READ(10,*) rec
  END DO

  !Lettura dei nodi della mesh
  CALL time2D_input_node(10)

  !Lettura delle incidenze degli elementi di bordo
  READ(10,*) rec
  CALL time2D_input_edges(10)

  !Chiusura del file mesh
  CLOSE(10)

  !Caratterizzazione degli elementi della mesh
  CALL time2D_elabora_input

  RETURN

END SUBROUTINE time2D_read_mesh
