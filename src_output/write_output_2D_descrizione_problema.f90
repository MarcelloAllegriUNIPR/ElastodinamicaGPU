SUBROUTINE write_output_2D_descrizione_problema(file_path)

  USE variable_2Dgeneral
  
  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABILI !!!!!!!!!!!!!!!!!!!!!!!!!!  

  !Input
  CHARACTER(len=200),INTENT(IN):: file_path

  !Variabili locali della subroutine
  INTEGER(kind=4):: ind_nodes, ind_element, grado

  CHARACTER(len=200):: type

  !!!!!!!!!!!!!!!!!!!!!!!!!! CORPO della SUBROUTINE !!!!!!!!!!!!!!!!!!

  OPEN(20000,FILE=TRIM(file_path))

  !Scrittura del tipo di contorno
  IF(list_elements(1)%nodes(1).EQ.&
       list_elements(number_elements)%nodes(2)) THEN
     WRITE(20000,*) '(c)',' CONTORNO CHIUSO'
  ELSE
     WRITE(20000,*) '(a)',' CONTORNO APERTO'
  END IF
  
  !Scrittura del tempo finale di analisi
  WRITE(20000,*) ' '
  WRITE(20000,*) T_fin,' TEMPO FINALE di ANALISI'
  
  !Scrittura del numero di intervalli temporali
  WRITE(20000,*) Nt,' NUMERO di INTERVALLI TEMPORALI'
  WRITE(20000,*) ' '

  !Scrittura dei parametri fisici
  WRITE(20000,*) rho, ' rho (DENSITA'' di MASSA)'
  WRITE(20000,*) mu, ' mu (MODULO di SHEAR)'
  !WRITE(20000,*) nu, ' RAPPORTO di POISSON'
  WRITE(20000,*) lambda, ' lambda'
  WRITE(20000,*) velC_S, ' velocità cs '
  WRITE(20000,*) velC_P, ' velocità cp '
  WRITE(20000,*) ' '
  
  !Scrittura dei nodi della mesh
  DO ind_nodes=1,number_nodes
     WRITE(20000,*) list_nodes(ind_nodes)%coordinates(1),&
          list_nodes(ind_nodes)%coordinates(2),&
          list_nodes(ind_nodes)%coordinates(3),&
          ind_nodes
  END DO
  WRITE(20000,*) ' '

  !Scrittura degli elementi della mesh
  DO ind_element=1,number_elements
     IF(list_elements(ind_element)%type.EQ.'u') THEN
        type='DIRICHLET'
     ELSE
        type='NEUMANN'
     END IF 
     IF(list_elements(ind_element)%type.EQ.'u') THEN
        grado=grado_q
     ELSE
        grado=grado_u
     END IF   
     WRITE(20000,*) list_elements(ind_element)%nodes(1),&
          list_elements(ind_element)%nodes(2),&
          TRIM(type),grado, list_elements(ind_element)%ind_BD
  END DO
  
  CLOSE(20000)

  RETURN
  
END SUBROUTINE write_output_2D_descrizione_problema

