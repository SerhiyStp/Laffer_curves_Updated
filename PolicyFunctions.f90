module PolicyFunctions

    use Model_Parameters

    implicit none
    ! Policy functions for active period of life
    real(8), dimension (:,:,:,:,:,:,:,:,:,:), allocatable :: v 
    real(8), dimension (:,:,:,:,:,:,:,:,:,:), allocatable, target :: ev_aux ! Pseudo-value fn for married proposed by Lars, that tracks separately continuation values for single men and women
    real(8), dimension (:,:,:,:,:,:,:,:,:,:), allocatable, target :: v_aux ! Pseudo-value fn for married proposed by Lars, that tracks separately continuation values for single men and women
    real(8), dimension (:,:,:,:,:,:,:,:,:,:), allocatable :: ev
    real(8), dimension (:,:,:,:,:,:,:,:,:,:), allocatable :: c 
    real(8), dimension (:,:,:,:,:,:,:,:,:,:), allocatable :: k 
    real(8), dimension (:,:,:,:,:,:,:,:,:,:), allocatable :: nm,nf 
    
    ! LFP regime specific functions for married
    real(8), dimension (:,:,:,:,:,:,:,:,:,:,:,:), allocatable :: c_lfp, k_lfp, nm_lfp, nf_lfp
    real(8), dimension (:,:,:,:,:,:,:,:,:,:,:,:,:), allocatable :: v_lfp
    real(8), dimension (:,:,:,:,:,:,:,:,:,:), allocatable :: lfpm, lfpf      

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
    
    ! LFP regime specific functions for single
    real(8), dimension (:,:,:,:,:,:,:,:), allocatable :: cs_lfp, ks_lfp, vs_lfp, ns_lfp  
    real(8), dimension (:,:,:,:,:,:,:), allocatable :: lfps      

    real(8), dimension (:), allocatable :: c_grid
    real(8), dimension (:), allocatable :: wage_grid
    real(8), dimension (:), allocatable :: k_grid
    real(8), dimension (:,:), allocatable, target :: exp_grid
    real(8), dimension (:,:,:), allocatable :: laborm,laborf
    real(8), dimension (:,:), allocatable :: labormwork,laborfwork
    real(8), dimension (:,:), allocatable :: laborsinglem,laborsinglef
    ! Policy functions for retired
    real(8), dimension (:,:,:,:,:,:), allocatable :: c_ret, v_ret
    real(8), dimension (:,:,:,:,:,:), allocatable, target :: ev_ret 
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

    ! Object-oriented policy functions
    integer, parameter :: nx = nk
    integer, parameter :: ny = nexp
    integer, parameter :: nz = nexp
    integer, parameter :: kxx = 3 !4 !3
    integer, parameter :: kyy = 3 !4 !3
    integer, parameter :: kzz = 3 !4 !3
    
    type policy_fn_2d
        real(8) :: coefs(nx, ny)
        integer :: inbvx
        integer :: inbvy
        integer :: iloy
        real(8) :: tx(nx+kxx)
        real(8) :: ty(ny+kyy)
    contains
        procedure :: eval => eval_2d
        procedure :: set => set_2d
    end type policy_fn_2d    
    
    type policy_fn_3d
        real(8) :: fvals(nx, ny, nz)
        real(8) :: coefs(nx, ny, nz)
        integer :: inbvx
        integer :: inbvy
        integer :: inbvz
        integer :: iloy
        integer :: iloz        
        real(8) :: tx(nx+kxx)
        real(8) :: ty(ny+kyy)
        real(8) :: tz(nz+kzz)
        real(8) :: xgrid(nx)
        real(8) :: ygrid(ny)
        real(8) :: zgrid(nz)
    contains
        procedure :: lin_interp => lin_interp_3d
        procedure :: set => set_3d
        procedure :: reset => reset_3d
        procedure :: eval => eval_3d
    end type policy_fn_3d


    type (policy_fn_3d), dimension(na, nu, na, nu, nfc, nfcm, 2), target :: pol_v_aux, pol_ev_aux
    !type (policy_fn_3d), dimension(na, nu, na, nu, T, nfc, nfcm, 2, 2, 2), target :: pol_v_mar_lfp
    type (policy_fn_3d), allocatable, target :: pol_v_mar_lfp(:, :, :, :, :, :, :, :, :, :)
    type (policy_fn_2d), dimension(2, na, nu, T, nfc, 2) :: pol_v_sing_lfp
    
    type (policy_fn_3d), dimension(na, na), target :: EV_mar_ret_pf
    
contains
    
    ! 2-d functions
    subroutine set_2d(self, fvals, xgrid, ygrid)
        use bspline_sub_module, only: db2ink
        class(policy_fn_2d) :: self
        real(8) :: xgrid(:), ygrid(:)
        real(8), intent(in) :: fvals(nx, ny)
        integer :: iflag
        integer :: iknot
        
        iknot = 0
        call db2ink(xgrid, nx, ygrid, ny, fvals, &
            kxx, kyy, iknot, self%tx, self%ty, self%coefs, iflag) 
    end subroutine set_2d
    
    function eval_2d(self, x)
        use bspline_sub_module, only: db2val
        class(policy_fn_2d) :: self
        real(8), intent(in) :: x(2)
        real(8) :: eval_2d
        integer :: iflag
        real(8) :: w1(kyy)
        real(8) :: w0(3*max(kxx,kyy) )
        integer :: inbvx, inbvy, iloy
        
        inbvx = 1
        inbvy = 1
        iloy = 1        
        call db2val(x(1), x(2), 0, 0, &
                    self%tx, self%ty, nx, ny, kxx, kyy, &
                    self%coefs, eval_2d, iflag, &
                    inbvx, inbvy, iloy, w1, w0, extrap=.true.)         
    end function eval_2d    

    ! 3-d functions
    subroutine set_3d(self, fvals, xgrid, ygrid, zgrid)
        use bspline_sub_module, only: db3ink
        class(policy_fn_3d) :: self
        real(8) :: xgrid(:), ygrid(:), zgrid(:)
        real(8), intent(in) :: fvals(nx, ny, nz)
        integer :: iflag, iknot
        
        iknot = 0
        call db3ink(xgrid, nx, ygrid, ny, zgrid, nz, &
                    fvals, kxx, kyy, kzz, iknot, self%tx, self%ty, self%tz, &
                    self%coefs, iflag)      
        call self%reset()
    end subroutine set_3d
    
    subroutine reset_3d(self)
        class(policy_fn_3d) :: self
        
        self%inbvx = 1
        self%inbvy = 1
        self%inbvz = 1
        self%iloy = 1
        self%iloz = 1         
    
    end subroutine reset_3d    

    function eval_3d(self, x)
        use bspline_sub_module, only: db3val
        use Utilities, only: trilin_interp
        class(policy_fn_3d) :: self
        real(8), intent(in) :: x(3)
        real(8) :: eval_3d
        integer :: iflag
        integer :: inbvx, inbvy, inbvz, iloy, iloz
        real(8) :: ww2(kyy,kzz), ww1(kzz), ww0(3*max(kxx,kyy,kzz))

        !inbvx = 1
        inbvx = self%inbvx
        !inbvy = 1
        inbvy = self%inbvy
        !inbvz = 1
        inbvz = self%inbvz
        !iloy = 1        
        iloy = self%iloy
        !iloz = 1       
        iloz = self%iloz       
        
        call db3val(x(1), x(2), x(3), 0, 0, 0, &
                    self%tx, self%ty, self%tz, &
                    nx, ny, nz, kxx, kyy, kzz, &
                    self%coefs, eval_3d, iflag,&
                    inbvx, inbvy, inbvz, iloy, iloz, &
                    ww2, ww1, ww0, extrap=.true.)         
        
        
    end function eval_3d
    
    
    function lin_interp_3d(self, x) result(f)
        use Utilities, only: trilin_interp
        class(policy_fn_3d) :: self
        real(8), intent(in) :: x(3)
        real(8) :: f
        
        f = trilin_interp(self%xgrid, self%ygrid, self%zgrid, self%fvals, nx, ny, nz, x)
    end function lin_interp_3d

end module PolicyFunctions
