subroutine Vuextra_gpu(li,lj,i,j,hk,indice_i,indice_j)
    !
    !     INTEGRALE DELLE FUNZIONI DI FORMA li,lj SUGLI ELEMENTI i,j
    !     ELEMENTI SEPARATI
    !
      USE variable_2Dgeneral
      use cudafor
      use kernels
      use VuExtraGpu
      IMPLICIT NONE
    
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABILI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
      !Input
      INTEGER(kind=4),INTENT(IN):: i, j, li, lj, hk, indice_i, indice_j
      double precision :: result    
      
      double precision, dimension(:,:), device :: CalculationResults(10,4)
      
      !!!!!!!!!!!!!!!!!!!!!!!!!! CORPO della SUBROUTINE !!!!!!!!!!!!!!!!!!
    
        CalculationResults=0.d0
        
        punto_m_1=list_nodes(j)%coordinates(1:2)
        punto_m_2=list_nodes(j+1)%coordinates(1:2)
        punto_m_tilde_1=list_nodes(i)%coordinates(1:2)
        punto_m_tilde_2=list_nodes(i+1)%coordinates(1:2)
        
        estremo_m=sqrt((punto_m_2(1)-punto_m_1(1))**2+(punto_m_2(2)-punto_m_1(2))**2)        
        estremo_m_tilde=sqrt((punto_m_tilde_2(1)-punto_m_tilde_1(1))**2+(punto_m_tilde_2(2)-punto_m_tilde_1(2))**2)
        CA=punto_m_tilde_1(1)-punto_m_1(1)
        CB=(punto_m_tilde_2(1)-punto_m_tilde_1(1))/estremo_m_tilde
        CC=(punto_m_2(1)-punto_m_1(1))/estremo_m
        CD=punto_m_tilde_1(2)-punto_m_1(2)
        CE=(punto_m_tilde_2(2)-punto_m_tilde_1(2))/estremo_m_tilde
        CF=(punto_m_2(2)-punto_m_1(2))/estremo_m
        
        CBCFCECC=CB*CF-CE*CC
        CFCACCCD=CF*CA-CC*CD
        CBCDCECA=CB*CD-CE*CA
        CBCCCECF=CB*CC+CE*CF
        CACCCFCD=CA*CC+CF*CD
        CACBCDCE=CA*CB+CD*CE

        if((dabs(CBCFCECC).le.1.d-15).and.(dabs(CFCACCCD).le.1.d-15)) then         
          flag_extra=1
        elseif((dabs(CBCFCECC).le.1.d-15).and.(dabs(CFCACCCD).gt.1.d-15)) then
          flag_extra=2    !	paralleli         
        else
          flag_extra=3    !   non allineati non paralleli          
        endif

        call setInstanceCommonData(CA, CB, CC, CD, CE, CF, &
        flag_extra, estremo_m, li,lj, estremo_m_tilde, &
        CBCFCECC, CFCACCCD, CBCCCECF, CACCCFCD, CACBCDCE, CBCDCECA)
        
        ! print *, 'VuExtra:',i,j
        !pause

        call Vuextra_sub_gpu(CalculationResults,1)
        call Vuextra_sub_gpu(CalculationResults,2)
        call Vuextra_sub_gpu(CalculationResults,3)
        call Vuextra_sub_gpu(CalculationResults,4)
        
        !pause
        call VuextraReduction<<<VuExtraGrid,VuExtraBlock>>>(CalculationResults,i,j)
        !result = Vu_device(i,j)
        !print *, result,i,j
        !deallocate(CalculationResults)
    RETURN
END