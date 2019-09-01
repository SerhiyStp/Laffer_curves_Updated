    module PolicyFunctions
    
        use Model_Parameters
        use pppack
        
        implicit none
        ! Policy functions for active period of life
        real(8), dimension (:,:,:,:,:,:,:,:), allocatable :: v,vdum
        real(8), dimension (:,:,:,:,:,:,:,:), allocatable :: ev
        real(8), dimension (:,:,:,:,:,:,:,:), allocatable :: c,cdum
        real(8), dimension (:,:,:,:,:,:,:,:), allocatable :: k,gkdum
        real(8), dimension (:,:,:,:,:,:,:,:), allocatable :: nm,nf,nmdum,nfdum
        real(8), dimension (:,:,:,:,:,:,:,:), allocatable :: ev_spln_coefs
        real(8), dimension (:,:,:,:,:,:,:,:), allocatable :: v_spln_coefs
        
        real(8), dimension (:,:,:,:,:,:,:), allocatable :: vs
        real(8), dimension (:,:,:,:,:,:,:), allocatable :: evs
        real(8), dimension (:,:,:,:,:,:,:), allocatable :: evm !Expected value function in the case of marriage next period
        real(8), dimension (:,:,:,:,:,:,:), allocatable :: cs
        real(8), dimension (:,:,:,:,:,:,:), allocatable :: ks
        real(8), dimension (:,:,:,:,:,:,:), allocatable :: ns
        real(8), dimension (:,:,:,:,:,:,:), allocatable :: evs_spln_coefs, vs_spln_coefs
        real(8), dimension (:,:,:,:,:,:,:), allocatable :: evm_spln_coefs !Coefficients to interpolate expected value function in the case of marriage next period
        real(8), dimension (:,:,:,:,:,:), allocatable :: fpartner,fpartnerdum,fpartnerdum2
        real(8), dimension (:,:,:,:), allocatable :: mpartner,mpartnerdum,mpartnerdum2
        real(8), dimension (:,:), allocatable :: ability_prob
        
        real(8), dimension (:), allocatable :: c_grid
        real(8), dimension (:), allocatable :: wage_grid
        real(8), dimension (:), allocatable :: k_grid
        real(8), dimension (:,:), allocatable :: exp_grid
        real(8), dimension (:,:), allocatable :: INTERP2D
        real(8), dimension (:), allocatable :: exp_grid_dum
        real(8), dimension (:), allocatable :: K_KNOT
        real(8), dimension (:,:), allocatable :: EXP_KNOT
        real(8), dimension (:,:,:), allocatable :: laborm,laborf
        real(8), dimension (:,:), allocatable :: labormwork
        real(8), dimension (:,:), allocatable :: laborsinglem,laborsinglef
        ! Policy functions for retired
        real(8), dimension (:,:,:), allocatable :: ev_spln_coefs_ret,evs_spln_coefs_ret
        real(8), dimension (:,:,:), allocatable :: p_ev_spln_coefs_ret,p_evs_spln_coefs_ret  ! for testing only
        real(8), dimension (:,:), allocatable :: c_ret, v_ret
        real(8), dimension (:,:), allocatable :: k_ret
        real(8), dimension (:,:), allocatable :: vs_ret
        real(8), dimension (:,:), allocatable :: cs_ret
        real(8), dimension (:,:), allocatable :: ks_ret
        real(8), dimension (:), allocatable :: break
        
    contains

    subroutine ppp_csint(xgrid, ydata, spln_coefs, n)
        integer :: n
        real(8) :: xgrid(n)
        real(8) :: ydata(n)
        real(8) :: spln_coefs(4,n)
        real(8) :: c(4,n)
        integer :: ibcbeg, ibcend
        
        ibcbeg = 0
        ibcend = 0
        c(1,:) = ydata
        call cubspl(xgrid, c, n, ibcbeg, ibcend)
        spln_coefs = c
    
        !call d_csint(k_grid, v_ret(:,Tret-it+1), BREAK, ev_spln_coefs_ret(:,:,Tret-it+1))
    
    end subroutine ppp_csint

    function ppp_csval(x, xgrid, n, spln_coefs)
        real(8) :: x
        integer :: n
        real(8) :: xgrid(n)
        real(8) :: spln_coefs(4,n)
        real(8) :: ppp_csval
        
        ppp_csval = ppvalu(xgrid, spln_coefs, n-1, 4, x, 0)
        
    end function ppp_csval
    
        
    end module PolicyFunctions
