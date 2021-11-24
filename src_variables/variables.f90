module Variables_DP

    implicit none

    double precision, device :: cs_d, cp_d, rho_d, delta_t_d
    double precision, device, allocatable, dimension(:) :: estremo_m, estremo_m_tilde, CB, CC, CE, CF
    double precision, device, allocatable, dimension(:) :: wint_d, xint_d, w_d, x_d 
    double precision, device, allocatable, dimension(:,:) :: punto_m_1, punto_m_2, punto_m_tilde_1, punto_m_tilde_2, CA, CD, CBCFCECC, CFCACCCD, &
                                                             CBCDCECA, CBCCCECF, CACCCFCD, CACBCDCE, Vu_d, matrix_d
    logical, device, allocatable, dimension(:) :: IsULeftAssI, IsULeftAssJ
    integer, device :: delta_kronecker_d(2,2), grado_q_d, DimVu_d, hk_d
    integer, device, allocatable, dimension(:,:) :: flag_extra, list_elements_d
end module Variables_DP