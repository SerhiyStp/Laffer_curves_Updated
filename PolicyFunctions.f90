 module PolicyFunctions
    
        use Model_Parameters
        
        implicit none
        ! Policy functions for active period of life
        real(8), dimension (:,:,:,:,:,:,:,:,:,:), allocatable :: v 
        real(8), dimension (:,:,:,:,:,:,:,:,:,:), allocatable :: ev
        real(8), dimension (:,:,:,:,:,:,:,:,:,:), allocatable :: c 
        real(8), dimension (:,:,:,:,:,:,:,:,:,:), allocatable :: k 
        real(8), dimension (:,:,:,:,:,:,:,:,:,:), allocatable :: nm,nf 
        
        real(8), dimension (:,:,:,:,:,:,:), allocatable :: vs
        real(8), dimension (:,:,:,:,:,:,:), allocatable :: evs
        real(8), dimension (:,:,:,:,:,:,:), allocatable :: evm !Expected value function in the case of marriage next period
        real(8), dimension (:,:,:,:,:,:,:), allocatable :: cs, edcs, Uprimes 
        real(8), dimension (:,:,:,:,:,:,:), allocatable :: ks
        real(8), dimension (:,:,:,:,:,:,:), allocatable :: ns
        real(8), dimension (:,:,:,:,:,:), allocatable :: fpartner,fpartnerdum,fpartnerdum2
        real(8), dimension (:,:,:,:,:,:), allocatable :: mpartner,mpartnerdum,mpartnerdum2
        real(8), dimension (:,:), allocatable :: ability_prob
        real(8), dimension (:,:,:), allocatable :: av_earnings
        
        real(8), dimension (:), allocatable :: c_grid
        real(8), dimension (:), allocatable :: wage_grid
        real(8), dimension (:), allocatable :: k_grid
        real(8), dimension (:,:), allocatable, target :: exp_grid
        real(8), dimension (:,:,:), allocatable :: laborm,laborf
        real(8), dimension (:,:), allocatable :: labormwork,laborfwork
        real(8), dimension (:,:), allocatable :: laborsinglem,laborsinglef
        ! Policy functions for retired
        real(8), dimension (:,:,:,:,:,:), allocatable :: c_ret, v_ret, ev_ret 
        real(8), dimension (:,:,:,:,:,:), allocatable :: edc_ret, Uprime_ret 
        real(8), dimension (:,:,:,:,:,:), allocatable :: k_ret
        real(8), dimension (:,:,:,:,:), allocatable :: vs_ret, evs_ret, Eulers_ret 
        real(8), dimension (:,:,:,:,:), allocatable :: cs_ret 
        real(8), dimension (:,:,:,:,:), allocatable :: edcs_ret, Uprimes_ret 
        real(8), dimension (:,:,:,:,:), allocatable :: ks_ret
        real(8), dimension (:), allocatable :: break
        
         ! New splines:        
        integer, parameter :: kx = 3
        integer, parameter :: ky = 3
        integer, parameter :: kz = 3
        real(8) :: tx(nk+kx)
        real(8) :: ty(nexp+ky,T+Tret)
        real(8) :: tz(nexp+kz,T+Tret)
        type bs_coefs_2d
            real(8) :: coefs(nk, nexp)
        end type bs_coefs_2d
        type bs_coefs_3d
            real(8) :: coefs(nk, nexp, nexp)
        end type bs_coefs_3d        
        
            
        type(bs_coefs_2d) :: evs_ret_bspl(2, na)
        type(bs_coefs_2d) :: edcs_ret_bspl(2, na)

        type(bs_coefs_2d) :: evs_bspl(2,na,nu,nfc)
        type(bs_coefs_2d), target :: vs_bspl(2,na,nu,nfc)
        type(bs_coefs_2d) :: evm_bspl(2,na,nu,nfc)
        
        type(bs_coefs_3d) :: ev_ret_bspl(na,na)
        type(bs_coefs_3d) :: edc_ret_bspl(na,na)

        type(bs_coefs_3d) :: ev_bspl(na,nu,na,nu,nfc,nfcm)
        type(bs_coefs_3d), target :: v_bspl(na,nu,na,nu,nfc,nfcm)       
        
    end module PolicyFunctions
