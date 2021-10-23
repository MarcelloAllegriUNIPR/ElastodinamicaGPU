MODULE variable_2Dgeneral
    
    use cudafor    
    !-- gpu --
    integer :: useGpu = 1

    !!! FILE !!!
    CHARACTER(len=200)::  dir, ext
    CHARACTER(len=8):: fmt
    
  !!!!!!!!!!!!!!!!!!!!!!!!!! PARAMETRI FISICI !!!!!!!!!!!!!!!!!!!!!!!!

  !Densità di massa rho
  REAL(kind=8):: rho
  !Modulo di shear mu
  REAL(kind=8):: mu
  !Rapporto di Poisson nu
  REAL(kind=8):: nu
  !Parametro di Lamé lambda
  REAL(kind=8):: lambda
  !Velocità delle onde P
  REAL(kind=8):: velC_P
  !Velocità delle onde S
  REAL(kind=8):: velC_S

  !!!!!!!!!!!!!!!! PARAMETRI DISCRETIZZAZIONE TEMPORALE !!!!!!!!!!!!!!

  !Tempo finale di analisi
  REAL(kind=8):: T_fin
  !Numero di intervalli temporali
  INTEGER(kind=4):: Nt
  !Passo di discretizzazione temporal
  REAL(kind=8)::deltat
  
  !!!!!!!!!!!!!!!!!!!!!!!! GRADO FUNZ di FORMA !!!!!!!!!!!!!!!!!!!!!!!

  !Grado delle funzioni di forma per approssimare lo spostamento
  INTEGER(kind=4):: grado_u
  !Grado delle funzioni di forma per approssimare la trazione
  INTEGER(kind=4):: grado_q
  
  !!!!!!!!!!!!!!!!!!!!!!!! VARIABILI del SISTEMA !!!!!!!!!!!!!!!!!!!!!!!
  
  INTEGER(kind=4):: DimVu
  !REAL(kind=8),DIMENSION(:),ALLOCATABLE:: u_bar
  !REAL(kind=8),DIMENSION(:),ALLOCATABLE:: rhs  
  !REAL(kind=8),DIMENSION(:,:),ALLOCATABLE:: Vu
  !REAL(kind=8),DIMENSION(:,:),ALLOCATABLE:: matrix
  
  !Variabile in cui viene memorizzata la soluzione del sistema
  REAL(kind=8),DIMENSION(:),ALLOCATABLE:: sol

  !!!!!!!!!!!!!!!!!!!!!!!!!! DESCRIZIONE MESH !!!!!!!!!!!!!!!!!!!!!!!!
  
  !Definizione tipo per gli elementi di bordo
  TYPE Selement
     INTEGER(kind=4),DIMENSION(2):: nodes !Estremi dell'elemento
     REAL(kind=8),DIMENSION(3):: normal   !Vettore normale
     REAL(kind=8):: length                !Lunghezza dell'elemento
     INTEGER(kind=4):: ind_BD             !Dato di bordo
     CHARACTER(len=1):: type              !Tipo di segmento di bordo
  END TYPE Selement

  TYPE(Selement),POINTER:: list_elements(:) => null()
  INTEGER(kind=4):: number_elements

  !Definizione tipo per i nodi della mesh
  TYPE Snode
     REAL(kind=8),DIMENSION(3):: coordinates !Coordinate dei nodi
     INTEGER(kind=4):: ind_domain            !Indice del sottodominio di
                                             !appartenenza
     INTEGER(kind=4):: eq                    !Indice per stabilire se il
                                             !nodo è portatore di incognita
  END TYPE Snode
  
  TYPE(Snode),POINTER:: list_nodes(:) => null()
  INTEGER(kind=4):: number_nodes
  
  !!!!!!!!!!!!!! FORMULE di QUADRATURA (GAUSS-LEGENDRE) !!!!!!!!!!!!!!
  
  TYPE gauss_type
      REAL(kind=8),DIMENSION(:),ALLOCATABLE:: nodiquad !Nodi di quadratura
      REAL(kind=8),DIMENSION(:),ALLOCATABLE:: pesiquad !Pesi di quadratura   
  END TYPE gauss_type
  
  TYPE(gauss_type),POINTER:: gauss(:) => null()
  INTEGER(kind=4):: number_gauss 

  INTEGER(kind=4) :: delta_kronecker_VuExtra(2,2), flag_extra, sharedMemDimension
  REAL(kind=8),DIMENSION(:),ALLOCATABLE :: x_VuExtra

  REAL(kind=8), DIMENSION(2) :: r, punto_m_1, punto_m_2, punto_m_tilde_1, punto_m_tilde_2
  double precision :: estremo_m, estremo_m_tilde, coeff_delta_kronecker, CA, CB, CC, CE, CF, CD,&
   CBCFCECC, CFCACCCD, CBCCCECF, CACCCFCD, CACBCDCE, CBCDCECA, cputime = 0.d0, gputime = 0.d0
END MODULE variable_2Dgeneral
